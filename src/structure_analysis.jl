using Graphs
using Random
using Plots
using Colors: RGB

"""
segments_from_backbone(x, y, edges)
Returns: segments, endpoints, junctions
"""
function segments_from_backbone(x::AbstractVector{<:Real},
                                y::AbstractVector{<:Real},
                                edges::Vector{Tuple{Int,Int}})
    N = length(x)
    @assert length(y) == N "x and y must have same length"

    g = SimpleGraph(N)
    for (i,j) in edges
        @assert 1 ≤ i ≤ N && 1 ≤ j ≤ N "edge endpoint out of bounds: ($i,$j)"
        i == j && continue             # ignore self-loops defensively
        add_edge!(g, i, j)
    end

    deg = degree.(Ref(g), 1:N)
    endpoints = [v for v in 1:N if deg[v] == 1]
    junctions = [v for v in 1:N if deg[v] ≥ 3]
    is_vertex = [deg[v] ≠ 2 for v in 1:N]   # endpoints or junctions

    # materialize neighbors for fast indexed access
    adj = [collect(neighbors(g, v)) for v in 1:N]

    # track undirected edges we’ve already traversed
    visited = Set{Tuple{Int,Int}}()
    und(u,v) = u < v ? (u,v) : (v,u)
    mark!(u,v) = push!(visited, und(u,v))
    seen(u,v) = und(u,v) in visited

    segments = Vector{Vector{Int}}()

    # Extend a path starting from directed edge u->v through a degree-2 chain.
    function extend_path(u::Int, v::Int; close_cycle::Bool=false)
        path = Int[u, v]
        mark!(u, v)
        first = u
        prev, cur = u, v
        while deg[cur] == 2
            n1, n2 = adj[cur][1], adj[cur][2]
            nxt = (n1 == prev) ? n2 : n1

            # if next edge already seen, stop (prevents infinite loop in cycles)
            if seen(cur, nxt)
                if close_cycle && nxt == first
                    push!(path, nxt)   # explicitly close the loop if requested
                end
                return path
            end
            push!(path, nxt)
            mark!(cur, nxt)
            prev, cur = cur, nxt
        end
        return path
    end

    # Traverse all branches starting from vertices (deg ≠ 2)
    for s in 1:N
        is_vertex[s] || continue
        for t in adj[s]
            seen(s, t) && continue
            push!(segments, extend_path(s, t))
        end
    end

    # Handle components that are pure cycles (no endpoints/junctions)
    if isempty(segments)
        for u in 1:N, v in adj[u]
            if !seen(u, v)
                loop = extend_path(u, v; close_cycle=true)
                push!(segments, loop)
            end
        end
    else
        # There can still be additional pure-cycle components in a multi-component graph
        for u in 1:N, v in adj[u]
            if !seen(u, v)
                loop = extend_path(u, v; close_cycle=true)
                push!(segments, loop)
            end
        end
    end

    return segments, endpoints, junctions
end


function segments_from_backbone_cc(x::AbstractVector, y::AbstractVector,
                                   edges::Vector{Tuple{Int,Int}})
    normalize_edges!(edges)
    g = SimpleGraph(max(length(x), length(y)))
    for (u,v) in edges; add_edge!(g, u, v); end
    segs = Vector{Vector{Int}}(); endpoints_all = Int[]; junctions_all = Int[]
    for comp in connected_components(g)
        setC = Set(comp)
        comp_edges = [(u,v) for (u,v) in edges if (u in setC && v in setC)]
        s,e,j = segments_from_backbone(x, y, comp_edges)  # your existing routine
        append!(segs, s); append!(endpoints_all, e); append!(junctions_all, j)
    end
    return segs, unique(endpoints_all), unique(junctions_all)
end


"""
plot_segments!(x, y, segments; endpoints=[], junctions=[])
Draw each segment in a distinct color; endpoints=cyan, junctions=green.
"""
function plot_segments!(x::AbstractVector, y::AbstractVector, segments::Vector{Vector{Int}};
                        endpoints::Vector{Int}=Int[], junctions::Vector{Int}=Int[])
    plt = plot(aspect_ratio=:equal, legend=false)

    # draw segments first
    rng = MersenneTwister(42)
    for seg in segments
        if length(seg) ≥ 2
            xs = x[seg]; ys = y[seg]
            plot!(xs, ys, lw=2, color=RGB(rand(rng), rand(rng), rand(rng)))
        end
    end

    # then draw vertex markers on top
    if !isempty(endpoints)
        scatter!(x[endpoints], y[endpoints], m=:circle, ms=6, color=:cyan)
    end
    if !isempty(junctions)
        scatter!(x[junctions], y[junctions], m=:circle, ms=6, color=:green)
    end
    return plt
end

############################
# Metrics helpers
############################

# Euclidean polyline length for a segment of vertex indices
function _polyline_length(x::AbstractVector, y::AbstractVector, seg::Vector{Int})
    L = length(seg)
    L ≤ 1 && return 0.0
    s = 0.0
    @inbounds for k in 1:L-1
        i, j = seg[k], seg[k+1]
        dx = x[j]-x[i]; dy = y[j]-y[i]
        s += hypot(dx,dy)
    end
    return s
end

# Per-segment curvature stats.
# Curvature = sum of absolute heading changes at interior vertices.
# Returns (total_abs_turn, mean_abs_turn_per_interior, curvature_density = total/length)
function segment_curvature_stats(x::AbstractVector, y::AbstractVector, seg::Vector{Int})
    L = length(seg)
    L < 3 && return (0.0, 0.0, 0.0)
    total = 0.0
    @inbounds for k in 2:L-1
        i0, i1, i2 = seg[k-1], seg[k], seg[k+1]
        v1x, v1y = x[i1]-x[i0], y[i1]-y[i0]
        v2x, v2y = x[i2]-x[i1], y[i2]-y[i1]
        n1 = hypot(v1x,v1y); n2 = hypot(v2x,v2y)
        (n1 == 0 || n2 == 0) && continue
        dot   = v1x*v2x + v1y*v2y
        cross = v1x*v2y - v1y*v2x
        # unsigned interior turning angle in [0, π]
        ang = atan(abs(cross), dot)
        total += ang
    end
    interior = max(L-2, 1)
    len = _polyline_length(x,y,seg)
    mean_abs = total / interior
    dens = len > 0 ? total / len : 0.0
    return (total, mean_abs, dens)
end

# Degree table from edges (undirected)
function degree_map(edges::Vector{Tuple{Int,Int}})
    deg = Dict{Int,Int}()
    @inbounds for (u,v) in edges
        deg[u] = get(deg,u,0) + 1
        deg[v] = get(deg,v,0) + 1
    end
    return deg
end

# Graph-based turning angle (deviation from straight) for deg==2 vertices
# Returns (indices::Vector{Int}, dev_angles_deg::Vector{Float64})
function graph_turn_deviation_deg(x::AbstractVector, y::AbstractVector,
                                  edges::Vector{Tuple{Int,Int}};
                                  min_edge_len::Real=1e-9)
    # adjacency
    adj = Dict{Int, Vector{Int}}()
    for (u,v) in edges
        push!(get!(adj,u,Int[]), v)
        push!(get!(adj,v,Int[]), u)
    end
    idxs = Int[]; devs = Float64[]
    for (i, nbrs) in adj
        length(nbrs) == 2 || continue
        a, b = nbrs[1], nbrs[2]
        v1x, v1y = x[a]-x[i], y[a]-y[i]
        v2x, v2y = x[b]-x[i], y[b]-y[i]
        n1 = hypot(v1x,v1y); n2 = hypot(v2x,v2y)
        (n1 < min_edge_len || n2 < min_edge_len) && continue
        ang = atan(abs(v1x*v2y - v1y*v2x), v1x*v2x + v1y*v2y)  # [0,π]
        dev = π - ang                                          # deviation from straight
        push!(idxs, i)
        push!(devs, dev * 180/π)
    end
    return idxs, devs
end

using StatsBase

# Plot and/or return segment length stats
function plot_segment_length_hist(x, y, segments; nbins::Int=30, save_path::Union{Nothing,String}=nothing)
    lens = [ _polyline_length(x,y,seg) for seg in segments ]
    plt = histogram(lens; bins=nbins, xlabel="segment length (plot units)", ylabel="freq",
                    title="Segment length distribution", legend=false)
    save_path !== nothing && savefig(plt, save_path)
    return lens, plt
end

# Plot turning-angle (graph-based deviation) histogram
function plot_turning_angle_hist(dev_angles_deg::Vector{Float64}; nbins::Int=30, save_path::Union{Nothing,String}=nothing)
    plt = histogram(dev_angles_deg; bins=nbins, xlabel="turning angle deviation (deg)", ylabel="freq",
                    title="Turning angles (graph-based)", legend=false)
    save_path !== nothing && savefig(plt, save_path)
    return plt
end

# Smooth density of turning angles (deviation-from-straight in degrees)
function plot_turning_angle_density(dev_angles_deg::Vector{Float64};
                                    npoints::Int=512,
                                    bandwidth::Union{Nothing,Real,Symbol}=nothing,
                                    save_path::Union{Nothing,String}=nothing)
    isempty(dev_angles_deg) && return nothing
    kd = isnothing(bandwidth) ? kde(dev_angles_deg) : kde(dev_angles_deg, bandwidth)
    xs, ys = kd.x, kd.density

    plt = plot(xs, ys; xlabel="turning angle deviation (deg)",
               ylabel="density", title="Turning angles (density)",
               legend=false, lw=2)
    scatter!(plt, dev_angles_deg, zero.(dev_angles_deg);
             markersize=2, markerstrokewidth=0, alpha=0.35, label=false)

    save_path !== nothing && savefig(plt, save_path)
    return plt
end


using KernelDensity   # Pkg.add("KernelDensity") if needed


# Example:
# x = state.x_coords; y = state.y_coords
# edges = collect(zip(1:length(x)-1, 2:length(x)))  # simple chain
# segments, endpoints, junctions = segments_from_backbone(x, y, edges)
# plt = plot_segments!(x, y, segments; endpoints=endpoints, junctions=junctions)
# savefig(plt, "backbone_segments.png")
