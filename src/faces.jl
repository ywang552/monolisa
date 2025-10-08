# === DCEL-lite face enumeration on a planar straight-line graph ===
# Works with your SimulationState.{x_coords, y_coords, edges}.
# Returns interior faces (vertex cycles) and their areas.

module Faces

export FaceResult, compute_enclosed_faces

struct FaceResult
    faces::Vector{Vector{Int}}      # each is a sequence of vertex ids (closed implicitly)
    areas::Vector{Float64}          # signed areas (CCW > 0)
    outer_face_idx::Int             # index in faces/areas of the dropped outer face (or 0)
    E_dedup::Vector{Tuple{Int,Int}} # deduplicated undirected edges actually used
end

# ------------------------- utilities -------------------------

@inline function poly_area_signed(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real})
    @assert length(xs) == length(ys)
    n = length(xs)
    acc = 0.0
    @inbounds for i in 1:n
        j = (i == n) ? 1 : (i+1)
        acc += xs[i]*ys[j] - xs[j]*ys[i]
    end
    return 0.5 * acc
end

# dedup and canonicalize undirected edges; remove self-loops
function dedup_edges(edges::Vector{Tuple{Int,Int}})
    seen = Set{Tuple{Int,Int}}()
    out  = Tuple{Int,Int}[]
    outreserve = sizehint!(out, length(edges))
    for (u,v) in edges
        if u == v; continue; end
        a,b = (u < v) ? (u,v) : (v,u)
        e = (a,b)
        if !(e in seen)
            push!(out, e); push!(seen, e)
        end
    end
    return out
end

# ------------------------- core: build half-edges -------------------------

"""
Build DCEL-lite arrays for directed half-edges.

Returns:
- he_src::Vector{Int}  # size = 2E
- he_dst::Vector{Int}
- he_rev::Vector{Int}  # index of reverse half-edge
- out_he::Vector{Vector{Int}}  # per-vertex outgoing half-edge ids, CCW-sorted
"""
function build_halfedges(N::Int, x::Vector{Float64}, y::Vector{Float64},
                         E::Vector{Tuple{Int,Int}})

    # Create directed half-edges
    he_src = Int[]; he_dst = Int[]
    sizehint!(he_src, 2length(E)); sizehint!(he_dst, 2length(E))
    for (u,v) in E
        push!(he_src, u); push!(he_dst, v)
        push!(he_src, v); push!(he_dst, u)
    end
    H = length(he_src)                # number of directed half-edges
    @assert H == 2*length(E)

    # Map (u,v) -> half-edge id
    hmap = Dict{Tuple{Int,Int},Int}()
    sizehint!(hmap, H)
    @inbounds for hid in 1:H
        hmap[(he_src[hid], he_dst[hid])] = hid
    end

    # Reverse indices
    he_rev = Vector{Int}(undef, H)
    @inbounds for hid in 1:H
        u = he_src[hid]; v = he_dst[hid]
        he_rev[hid] = hmap[(v,u)]
    end

    # Angles per half-edge (for sorting outgoing lists)
    theta = Vector{Float64}(undef, H)
    @inbounds for hid in 1:H
        u = he_src[hid]; v = he_dst[hid]
        theta[hid] = atan(y[v]-y[u], x[v]-x[u])  # atan2(y, x)
    end

    # Build per-vertex outgoing half-edge lists, sorted CCW by theta
    out_he = [Int[] for _ in 1:N]
    @inbounds for hid in 1:H
        push!(out_he[he_src[hid]], hid)
    end
    @inbounds for u in 1:N
        sort!(out_he[u], by = h -> theta[h])    # CCW order
    end

    return he_src, he_dst, he_rev, out_he
end

# ------------------------- face walking -------------------------

"""
Given a directed half-edge id h0, walk the face by always taking the "next-left" turn.
Requires:
- he_src, he_dst, he_rev
- out_he: CCW-sorted outgoing half-edges per vertex
- pos_in_out: cache of positions of each half-edge in its source vertex's out_he list (optional speedup)

Returns the sequence of vertex ids forming the boundary (closed implicitly),
and the list of half-edges visited in this loop.
"""
function walk_face(h0::Int,
                   he_src::Vector{Int}, he_dst::Vector{Int},
                   he_rev::Vector{Int}, out_he::Vector{Vector{Int}},
                   pos_in_out::Union{Nothing,Vector{Int}}=nothing)

    h = h0
    verts = Int[]         # boundary vertices
    hedges = Int[]        # directed half-edges in this loop

    while true
        push!(hedges, h)
        u = he_src[h]; v = he_dst[h]
        push!(verts, u)

        # at vertex v, find reverse half-edge index and take next CCW neighbor
        rev = he_rev[h]          # this is v->u
        lst = out_he[v]
        # locate rev in lst:
        if pos_in_out === nothing
            idx = findfirst(==(rev), lst)
        else
            idx = pos_in_out[rev]  # O(1)
        end
        @assert idx !== nothing
        i = Int(idx)
        # i_next = (i == length(lst)) ? 1 : i+1

        i_next = (i == 1) ? length(lst) : i-1
        h = lst[i_next]           # next-left half-edge: v -> w

        if h == h0
            # closed the loop; add final vertex if you want explicit closure (u_1 again)
            break
        end
    end
    return verts, hedges
end

# Precompute positions of each half-edge in its source-vertex outgoing list
function index_positions(out_he::Vector{Vector{Int}})
    Hmax = 0
    for lst in out_he
        if !isempty(lst)
            Hmax = max(Hmax, maximum(lst))
        end
    end
    pos = fill(0, Hmax)
    for u in 1:length(out_he)
        lst = out_he[u]
        @inbounds for (k,h) in enumerate(lst)
            pos[h] = k
        end
    end
    return pos
end

# ------------------------- main API -------------------------

"""
compute_enclosed_faces(state; area_floor=0.0, drop_outer=true)

- Enumerates all faces via angular-walk.
- Computes signed areas (CCW positive).
- Drops outer face if requested (largest |area|).
- Filters faces with |area| < area_floor.

Returns FaceResult.
"""
# --- replace inside compute_enclosed_faces(...) ---

function compute_enclosed_faces(state;
        area_floor::Float64 = 0.0,
        drop_outer::Bool = true,
        normalize_orientation::Bool = true,
        return_abs::Bool = true)

    x = state.x_coords
    y = state.y_coords
    N = min(length(x), length(y))
    @assert N > 0 "No placed monomers."
    if length(x) != length(y)
        @warn "x/y length mismatch: x=$(length(x)), y=$(length(y)). Truncating to N=$N."
        x = view(x, 1:N); y = view(y, 1:N)
    end

    ## TODO delete this later
    rawE = dedup_edges(state.edges)
    E = Tuple{Int,Int}[]
    sizehint!(E, length(rawE))
    dropped = 0
    for (u,v) in rawE
        if 1 <= u <= N && 1 <= v <= N
            push!(E, (u,v))
        else
            dropped += 1
        end
    end
    dropped > 0 && @warn "Dropped $dropped edges referencing > N=$N."
    ####

    he_src, he_dst, he_rev, out_he = build_halfedges(N, x, y, E)
    pos = index_positions(out_he)

    H = length(he_src)
    visited = falses(H)
    faces = Vector{Vector{Int}}()
    areas = Float64[]

    @inbounds for h0 in 1:H
        if visited[h0]; continue; end
        verts, hedges = walk_face(h0, he_src, he_dst, he_rev, out_he, pos)
        for h in hedges; visited[h] = true; end
        xs = @view x[verts]; ys = @view y[verts]
        A = poly_area_signed(xs, ys)      # signed area; CCW > 0
        push!(faces, verts); push!(areas, A)
    end

    # Identify the outer face BEFORE orientation normalization / abs
    outer_idx = 0
    if !isempty(areas)
        _, outer_idx = findmax(abs.(areas))   # largest magnitude area is outer
    end

    # Optionally normalize orientation: make interior loops CCW
    if normalize_orientation
        for i in eachindex(faces)
            if areas[i] < 0.0
                reverse!(faces[i])            # flip vertex order
                areas[i] = -areas[i]          # make positive
            end
        end
    end

    # Drop outer and tiny slivers
    keep = trues(length(faces))
    if drop_outer && outer_idx != 0
        keep[outer_idx] = false
    end
    if area_floor > 0
        for i in eachindex(keep)
            keep[i] && abs(areas[i]) < area_floor && (keep[i] = false)
        end
    end

    faces_keep = [faces[i] for i in eachindex(faces) if keep[i]]
    areas_keep = [areas[i] for i in eachindex(areas) if keep[i]]

    # Optionally return absolute values for convenience
    if return_abs
        for i in eachindex(areas_keep)
            areas_keep[i] = abs(areas_keep[i])
        end
    end

    return FaceResult(faces_keep, areas_keep, drop_outer ? outer_idx : 0, E)
end


end # module Faces

using .Faces

# Compute CDFs from a Faces.FaceResult
struct FaceCDF
    x::Vector{Float64}   # sorted area thresholds
    y_count::Vector{Float64}   # Count-CDF: (#faces ≤ x) / (#faces)
    y_area::Vector{Float64}    # Area-share CDF: (sum areas ≤ x) / (total area)
end

function face_cdfs(res; min_area::Float64=0.0)
    a = abs.(res.areas)
    keep = findall(ai -> ai ≥ min_area, a)
    isempty(keep) && return FaceCDF(Float64[], Float64[], Float64[])
    a = sort(a[keep])  # ascending
    n = length(a)
    y_count = collect(1:n) ./ n
    s = cumsum(a)
    y_area = s ./ s[end]
    return FaceCDF(a, y_count, y_area)
end

"""
plot_face_area_cdf(res; logx=false, min_area=0.0, title="Face-size CDFs")

- Plots two curves on the same axes:
  * Count-CDF (fraction of faces)
  * Area-share CDF (fraction of total enclosed area)
- y ∈ [0,1]; x = face area (optionally log-scaled).
"""
function plot_face_area_cdf(res; logx::Bool=false, min_area::Float64=0.0,
                            title::AbstractString="Face-size CDFs", label = "C2sep", p = nothing)
    cdf = face_cdfs(res; min_area=min_area)
    if isempty(cdf.x)
        @warn "No faces after filtering (min_area=$(min_area))."
        return plot()
    end

    if isnothing(p)
        plt = plot(; legend=:bottomright, xlabel="Face area (nm²)",
                   ylabel="Cumulative fraction", ylims=(0,1),
                   aspect_ratio=:auto, title=title, grid=true)
    else
        plt = p

    end 

    if logx
        # Avoid log(0): shift the minimum > 0 slightly if needed
        xmin = maximum([minimum(cdf.x), nextfloat(0.0)])
        plot!(plt, cdf.x, cdf.y_count; xscale=:log10, label="Count-CDF $(label)", lw=2)
        # plot!(plt, cdf.x, cdf.y_area;  xscale=:log10, label="Area-share CDF", lw=2, ls=:dash)
        xlims!(1e0, 1e5)
    else
        plot!(plt, cdf.x, cdf.y_count; label="Count-CDF", lw=2)
        plot!(plt, cdf.x, cdf.y_area;  label="Area-share CDF", lw=2, ls=:dash)
    end
    return plt
end

# # res = Faces.compute_enclosed_faces(state; drop_outer=true, normalize_orientation=true, return_abs=true)
# res = Faces.compute_enclosed_faces(state; area_floor=5., drop_outer=true)



# plt_cdf = plot_face_area_cdf(res; logx=true, min_area=5.0,
#     title="Face-size CDFs")

# res2 = Faces.compute_enclosed_faces(state2; area_floor=5., drop_outer=true)
# plot_face_area_cdf(res2; logx = true, p = plt_cdf, min_area=5.0, label = "8C")

# savefig("plots\\tmp\\faceCDF.png")
