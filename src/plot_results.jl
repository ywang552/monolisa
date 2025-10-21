############################
# Monomer Plotting (units + spin tick)
############################



# -------------------------
# Config: physical mapping
# -------------------------
const NM_PER_DATA = 1.55 

# -------------------------
# Small utilities
# -------------------------
ceil_to(x, step) = step * ceil(Int, x / step)
build_ticks(limit, step) = collect(0:step:limit)
const DEG2RAD = π/180



"""
draw_edges!(plt, x, y, edges; lc=:grey, lw=2, alpha=1.0)
Overlay straight segments for each (i,j) pair.
"""
function draw_edges!(plt, x::AbstractVector, y::AbstractVector,
                     edges::Vector{Tuple{Int,Int}}; lc=:grey, lw=2, alpha=1.0)
    @inbounds for (i,j) in edges
        plot!(plt, @view(x[[i,j]]), @view(y[[i,j]]);
              seriestype=:path, color=lc, linewidth=lw, alpha=alpha)
    end
    return plt
end


# --- f1: features from x,y,edges ---------------------------------------------

"""
backbone_features(x, y, edges) -> (; segments, endpoints, junctions)

Thin wrapper around `segments_from_backbone` that returns a NamedTuple,
so the result can be passed around cleanly.
"""
function backbone_features(x::AbstractVector{<:Real},
                           y::AbstractVector{<:Real},
                           edges::Vector{Tuple{Int,Int}})
    segments, endpoints, junctions = segments_from_backbone(x, y, edges)
    return (; segments, endpoints, junctions)
end

# --- f2: draw the raw backbone edges -----------------------------------------

"""
draw_backbone!(plt, x, y, edges; lc=:grey, lw=2, alpha=1.0)

Draw straight segments for each undirected edge (i,j).
Does not compute features or color by segment—kept minimal by design.
"""
function draw_backbone!(plt,
                        x::AbstractVector, y::AbstractVector,
                        edges::Vector{Tuple{Int,Int}}; lc=:grey, lw::Real=2, alpha::Real=1.0)
    @inbounds for (i, j) in edges
        # discontiguous → make a tiny container; no @view here
        plot!(plt, [x[i], x[j]], [y[i], y[j]];
              seriestype=:path, color=lc, linewidth=lw, alpha=alpha)
    end
    return plt
end

# --- f3: overlay segments + vertex markers -----------------------------------

"""
overlay_segments_vertices!(plt, x, y, features;
                           segment_mode=:random, endpoints_color=:cyan,
                           junctions_color=:green, seg_lw=2)

Overlay per-segment polylines (distinct colors) and vertex markers (endpoints/junctions).
`features` is the NamedTuple returned by `backbone_features`.
segment_mode = :random (RGB via seeded RNG) or :cmap (categorical colormap).
"""
function overlay_segments_vertices!(plt,
                                    x::AbstractVector, y::AbstractVector,
                                    features::NamedTuple;
                                    segment_mode::Symbol = :random,
                                    endpoints_color=:cyan,
                                    junctions_color=:green,
                                    seg_lw::Real=2,
                                    v_ms = 2, 
                                    # NEW:
                                    turning_vertices::AbstractVector{<:Integer}=Int[],
                                    turning_color=:orange,
                                    turning_ms::Real=6,
                                    turning_marker::Symbol=:square)

    segments    = features.segments
    endpoints   = get(features, :endpoints, Int[])
    junctions   = get(features, :junctions, Int[])

    if segment_mode === :random
        rng = MersenneTwister(42)
        for seg in segments
            if length(seg) ≥ 2
                plot!(plt, x[seg], y[seg];
                      lw=seg_lw, color=RGB(rand(rng), rand(rng), rand(rng)))
            end
        end
    elseif segment_mode === :cmap
        cmap = cgrad(:viridis, max(length(segments), 1), categorical=true)
        for (k, seg) in enumerate(segments)
            if length(seg) ≥ 2
                plot!(plt, x[seg], y[seg]; lw=seg_lw, color=cmap[k])
            end
        end
    else
        error("segment_mode must be :random or :cmap")
    end

    # draw markers last so they sit on top
    if !isempty(endpoints)
        scatter!(plt, x[endpoints], y[endpoints]; m=:circle, ms=v_ms, color=endpoints_color)
    end
    if !isempty(junctions)
        scatter!(plt, x[junctions], y[junctions]; m=:circle, ms=v_ms, color=junctions_color)
    end
    if !isempty(turning_vertices)
        scatter!(plt, x[turning_vertices], y[turning_vertices];
                 m=:square, ms=turning_ms, color=turning_color)
    end
    return plt
end

"""
Compare brute-force edges (from find_contact_edges) with direct edges (state.edges).
Prints counts and mismatches.
"""
function compare_edges(state, edges_bruteforce::Vector{Tuple{Int,Int}})
    # normalize to undirected, no self-loops
    norm(edges) = Set([(min(u,v), max(u,v)) for (u,v) in edges if u != v])

    E_brute  = norm(edges_bruteforce)
    E_direct = norm(state.edges)

    missing_in_direct = setdiff(E_brute, E_direct)
    extra_in_direct   = setdiff(E_direct, E_brute)

    println("Brute-force edges: ", length(E_brute))
    println("Direct edges:      ", length(E_direct))
    println("Missing in direct: ", length(missing_in_direct))
    println("Extra in direct:   ", length(extra_in_direct))
    for (u,v) in Iterators.take(missing_in_direct, 5)
        dx = state.x_coords[u] - state.x_coords[v]
        dy = state.y_coords[u] - state.y_coords[v]
        d  = hypot(dx,dy)
        println((u,v), "  d=", d, "  2r=", 2*state.radius, "  2r*1.02=", 2*state.radius*CONTACT_SCALE)
    end

    return (; E_brute, E_direct, missing_in_direct, extra_in_direct)
end


function _draw_contacts_grid!(plt, x::AbstractVector, y::AbstractVector, r::Real;
                              xlim::Tuple, ylim::Tuple,
                              contact_scale::Real=1.02,
                              max_lines::Int=200_000, lw=2, color=:black)
    edges = find_contact_edges(x, y, r; xlim=xlim, ylim=ylim,
                               contact_scale=contact_scale, max_edges=max_lines)
    draw_backbone!(plt, x, y, edges; lc=color, lw=lw, alpha=1.0)
    return edges
end


function find_contact_edges(x::AbstractVector, y::AbstractVector, r::Real;
                            xlim::Tuple, ylim::Tuple, contact_scale::Real=1.02,
                            max_edges::Int=200_000)
    @assert length(x) == length(y)
    N = length(x)
    contact_r = 2 * r * contact_scale   # <-- fix 2*r
    csq = contact_r^2
    cell = contact_r

    gx(i) = Int(floor((x[i] - first(xlim)) / cell))
    gy(i) = Int(floor((y[i] - first(ylim)) / cell))

    buckets = Dict{Tuple{Int,Int}, Vector{Int}}()
    @inbounds for i in 1:N
        xi, yi = x[i], y[i]
        (isfinite(xi) && isfinite(yi)) || continue
        push!(get!(buckets, (gx(i), gy(i)), Int[]), i)
    end

    offs  = ((-1,-1),(-1,0),(-1,1),(0,-1),(0,0),(0,1),(1,-1),(1,0),(1,1))
    edges = Vector{Tuple{Int,Int}}()

    @inbounds for (gi, gj) in keys(buckets)
        ids = buckets[(gi,gj)]
        for (di, dj) in offs
            nbrs = get(buckets, (gi+di, gj+dj), nothing)
            isnothing(nbrs) && continue
            for i in ids
                xi, yi = x[i], y[i]
                for j in nbrs
                    j <= i && continue
                    dx = x[j]-xi; dy = y[j]-yi
                    if dx*dx + dy*dy <= csq
                        push!(edges, (i,j))
                        length(edges) >= max_edges && return edges
                    end
                end
            end
        end
    end
    return edges
end

# Minimal Union-Find for Kruskal MST
mutable struct UF
    parent::Vector{Int}
    rank::Vector{Int}
end
UF(n::Int) = UF(collect(1:n), fill(0, n))
function findset(uf::UF, a::Int)
    while uf.parent[a] != a
        uf.parent[a] = uf.parent[uf.parent[a]]
        a = uf.parent[a]
    end
    return a
end
function unite!(uf::UF, a::Int, b::Int)
    ra, rb = findset(uf, a), findset(uf, b)
    ra == rb && return false
    if uf.rank[ra] < uf.rank[rb]
        uf.parent[ra] = rb
    elseif uf.rank[rb] < uf.rank[ra]
        uf.parent[rb] = ra
    else
        uf.parent[rb] = ra
        uf.rank[ra] += 1
    end
    return true
end

"""
    backbone_edges_mst(x, y, r; xlim, ylim, contact_scale=1.02)

Return a *minimum-spanning forest* (MST per connected component) over the
contact graph. Only returns the edges (i, j) — no plotting.
"""
function backbone_edges_mst(x::AbstractVector, y::AbstractVector, r::Real;
                            xlim::Tuple, ylim::Tuple, contact_scale::Real=1.02)
    N = length(x)
    edges, dists = find_contact_edges(x, y, r; xlim=xlim, ylim=ylim, contact_scale=contact_scale)

    # Sort edges by distance
    order = sortperm(dists)
    uf = UF(N)
    backbone = Vector{Tuple{Int,Int}}()
    for k in order
        i, j = edges[k]
        unite!(uf, i, j) && push!(backbone, (i, j))
    end
    return backbone
end

"""
turning_vertices_from_edges(x, y, edges; angle_thresh_deg=30.0)

Mark nodes with degree==2 whose angle between the two incident edges
(at the node) is >= angle_thresh_deg. Returns Vector{Int}.
"""
function turning_vertices_from_edges(x::AbstractVector, y::AbstractVector,
                                     edges::Vector{Tuple{Int,Int}};
                                     angle_thresh_deg::Real=20.0,
                                     min_edge_len::Real=1e-9)

    th = deg2rad(angle_thresh_deg)

    # adjacency
    adj = Dict{Int, Vector{Int}}()
    for (u,v) in edges
        push!(get!(adj,u,Int[]), v)
        push!(get!(adj,v,Int[]), u)
    end

    turns = Int[]
    for (i, nbrs) in adj
        length(nbrs) == 2 || continue
        a, b = nbrs[1], nbrs[2]

        v1x, v1y = x[a]-x[i], y[a]-y[i]
        v2x, v2y = x[b]-x[i], y[b]-y[i]
        n1 = hypot(v1x, v1y); n2 = hypot(v2x, v2y)
        (n1 < min_edge_len || n2 < min_edge_len) && continue

        # angle between rays ∈ [0, π]
        dot   = v1x*v2x + v1y*v2y
        cross = v1x*v2y - v1y*v2x
        ang = atan(abs(cross), dot)      # 0..π

        # deviation from straight (π) → 0 for straight, large for turns
        dev = π - ang

        if dev ≥ th
            push!(turns, i)
        end
    end
    return turns
end


"""
turning_vertices_from_segments(x, y, segments;
                               angle_thresh_deg=60.0, min_gap=2)
Return a Set{Int} of nodes (degree-2 along the path) whose signed turn angle
|Δθ| exceeds angle_thresh_deg. min_gap prevents two splits too close together.
"""
function turning_vertices_from_segments(x::AbstractVector, y::AbstractVector,
                                        segments::Vector{Vector{Int}};
                                        angle_thresh_deg::Real=60.0,
                                        min_gap::Int=2)
    th = deg2rad(angle_thresh_deg)
    turning = Set{Int}()

    for seg in segments
        L = length(seg)
        L < 3 && continue
        last_split_k = -typemax(Int)
        @inbounds for k in 2:L-1
            i_prev, i, i_next = seg[k-1], seg[k], seg[k+1]
            v1x, v1y = x[i] - x[i_prev], y[i] - y[i_prev]
            v2x, v2y = x[i_next] - x[i], y[i_next] - y[i]
            # skip degenerate steps
            norm1 = hypot(v1x, v1y); norm2 = hypot(v2x, v2y)
            (norm1 == 0 || norm2 == 0) && continue

            # signed turn angle in (-π, π]
            a1 = atan(v1y, v1x)
            a2 = atan(v2y, v2x)
            dθ = a2 - a1
            dθ = (dθ + π) % (2π) - π  # wrap to (-π,π]

            if abs(dθ) ≥ th && (k - last_split_k) ≥ min_gap
                push!(turning, i)
                last_split_k = k
            end
        end
    end
    return turning
end

"""
refine_segments_by_curvature(x, y, segments;
    angle_thresh_deg=60.0, min_gap=2,
    keep_short=true, min_len=2,
    include_turn_in_prev=true)

- Splits at curvature peaks (from turning_vertices_from_segments).
- Ensures segments have >=2 nodes unless keep_short=true.
- include_turn_in_prev: if true, left piece ends at the turning vertex;
  if false, left piece ends at (k-1) and the turning vertex starts the right piece.
"""
function refine_segments_by_curvature(x::AbstractVector, y::AbstractVector,
                                      segments::Vector{Vector{Int}};
                                      angle_thresh_deg::Real=60.0,
                                      min_gap::Int=2,
                                      keep_short::Bool=true,
                                      min_len::Int=2,
                                      include_turn_in_prev::Bool=true)

    # detect turning nodes (degree-2 along path) above threshold
    tset = turning_vertices_from_segments(x, y, segments;
                                          angle_thresh_deg=angle_thresh_deg,
                                          min_gap=min_gap)

    new_segments = Vector{Vector{Int}}()

    # helper: push only if useful length
    @inline function maybe_push!(out, segslice)
        L = length(segslice)
        if L ≥ max(2, min_len) || (keep_short && L ≥ 2)
            push!(out, segslice)
        end
    end

    for seg in segments
        L = length(seg)
        if L < 3
            maybe_push!(new_segments, seg)
            continue
        end

        start_idx = 1
        k = 2
        while k ≤ L-1
            if seg[k] in tset
                # choose split boundary behavior
                left_end = include_turn_in_prev ? k : (k-1)
                right_start = k

                # left piece [start_idx .. left_end]
                if left_end ≥ start_idx
                    maybe_push!(new_segments, seg[start_idx:left_end])
                end

                # next piece starts at the turning vertex (or its index)
                start_idx = right_start
            end
            k += 1
        end

        # tail
        if start_idx ≤ L
            maybe_push!(new_segments, seg[start_idx:L])
        end
    end

    return new_segments, collect(tset)
end



# -------------------------
# Transform (shift + optional fit-scale)
# -------------------------
function transform_coords(x::AbstractVector, y::AbstractVector;
                          mode::Symbol=:min,
                          scale::Union{Symbol,Real}=:identity,
                          target_max::Union{Nothing,Real}=nothing)
    @assert length(x) == length(y) "x and y must have same length"

    # filter non-finite
    mask = isfinite.(x) .& isfinite.(y)
    x0 = collect(@view x[mask]); y0 = collect(@view y[mask])

    # handle empty safely
    if isempty(x0)
        return similar(x, 0), similar(y, 0), 1.0
    end

    # shift
    if mode === :min
        x0 .-= minimum(x0);  y0 .-= minimum(y0)
    elseif mode === :center
        x0 .-= (minimum(x0) + maximum(x0))/2
        y0 .-= (minimum(y0) + maximum(y0))/2
    else
        error("mode must be :min or :center")
    end

    # scale
    s = if scale === :identity
        1.0
    elseif scale === :fit
        isnothing(target_max) && error("need target_max for scale=:fit")
        bigx = maximum(x0) - minimum(x0)
        bigy = maximum(y0) - minimum(y0)
        big  = max(bigx, bigy)
        big <= 0 ? 1.0 : Float64(target_max)/big
    else
        Float64(scale)
    end

    # tiny padding so origin isn't at 0,0 corner
    x0 .+= 5
    y0 .+= 5

    return x0 .* s, y0 .* s, s
end

# -------------------------
# Unit helpers (auto-pick nm / μm / mm)
# -------------------------
pick_unit(span_nm::Real) = span_nm ≥ 1e6 ? :mm : (span_nm ≥ 1e3 ? :μm : :nm)
unit_divisor(u::Symbol) = u === :nm ? 1.0 : (u === :μm ? 1e3 : 1e6)  # nm → unit
unit_label(u::Symbol)   = u === :μm ? "μm" : String(u)

# “nice” step targeting ~8 ticks
function nice_tick_step(span::Real; target_ticks::Int=8)
    raw = span / max(target_ticks, 1)
    if raw <= 0
        return 1.0
    end
    exp10 = floor(log10(raw))
    base  = raw / 10.0^exp10
    m = base ≤ 1.0 ? 1.0 : base ≤ 2.0 ? 2.0 : base ≤ 5.0 ? 5.0 : 10.0
    return m * 10.0^exp10
end

# -------------------------
# Density grid (for massive LOD)
# -------------------------
function density_grid(x::AbstractVector, y::AbstractVector; bins=(400, 400))
    hx = fit(Histogram, x, bins=bins[1])
    hy = fit(Histogram, y, bins=bins[2])
    counts = zeros(length(hx.edges[1])-1, length(hy.edges[1])-1)
    for (xi, yi) in zip(x, y)
        ix = searchsortedfirst(hx.edges[1], xi) - 1
        iy = searchsortedfirst(hy.edges[1], yi) - 1
        (1 ≤ ix ≤ size(counts,1) && 1 ≤ iy ≤ size(counts,2)) && (counts[ix, iy] += 1)
    end
    return counts', hx.edges[1], hy.edges[1]  # transpose for heatmap orientation
end

# -------------------------
# Spin tick: ring segment
# -------------------------

"""
Overlay rim-segment spin markers on visible points.
- indices: absolute indices that were actually scattered (e.g. collect(idx))
"""
function add_orientation_arcs!(plt, x, y, rotation, r;
                               indices, every=1,
                               arc_span_deg=25.0,
                               thickness_frac=0.45,
                               color=:blue, alpha=0.95)
    invlog10(x) = 1 - log10(x) / 6   # since log10(1_000_000) = 6

    # stroke width in data units (fraction of radius)
    lw = 5 * invlog10(length(x))  # data-units linewidth

    nseg = 24  # resolution for smooth arc
    for k in 1:every:length(indices)
        i = indices[k]
        θ0 = (rotation[i] - arc_span_deg/2) * π/180
        θ1 = (rotation[i] + arc_span_deg/2) * π/180
        θs = range(θ0, θ1; length=nseg)
        xs = x[i] .+ r .* cos.(θs)
        ys = y[i] .+ r .* sin.(θs)

        plot!(plt, xs, ys;
              seriestype=:path,
              color=color, alpha=alpha,
              linewidth=lw,                     # stroke thickness
              linecap=:round)                   # prettier ends
    end
end


"""
Choose the least-occupied corner for an inset compass.

Returns: (:topleft | :topright | :bottomleft | :bottomright, cx, cy, r)
- cx, cy = center of the compass
- r      = compass radius (in plot units)
"""
function _choose_compass_corner(xT::AbstractVector, yT::AbstractVector,
                                xlim::Tuple, ylim::Tuple;
                                size_frac::Real=0.12,    # compass radius relative to min span
                                margin_frac::Real=0.05,  # gap to axes
                                sample_cap::Int=100_000)

    N = length(xT)
    # sample for speed if huge
    if N > sample_cap
        stride = cld(N, sample_cap)
        xS = @view xT[1:stride:end]; yS = @view yT[1:stride:end]
    else
        xS = xT; yS = yT
    end

    xlo, xhi = xlim; ylo, yhi = ylim
    spanx = xhi - xlo; spany = yhi - ylo
    spanu = min(spanx, spany)

    r = size_frac * spanu
    m = margin_frac * spanu
    box = 2r   # square that bounds the compass

    corners = Dict(
        :topleft     => (xlo + m,         yhi - m - box),
        :topright    => (xhi - m - box,   yhi - m - box),
        :bottomleft  => (xlo + m,         ylo + m),
        :bottomright => (xhi - m - box,   ylo + m),
    )

    # Count points inside each corner box quickly
    counts = Dict{Symbol,Int}()
    for (k, (xb, yb)) in corners
        xc1, xc2 = xb, xb + box
        yc1, yc2 = yb, yb + box
        # boolean mask count (vectorized & fast)
        c = count(((xS .>= xc1) .& (xS .<= xc2) .& (yS .>= yc1) .& (yS .<= yc2)))
        counts[k] = c
    end

    # take the corner with the fewest points; in a tie prefer top-right
    order = sort(collect(keys(counts)); by=k->(counts[k], k==:topright ? 0 : 1))
    best = first(order)
    xb, yb = corners[best]

    # center of the compass circle for plotting
    cx, cy = xb + r, yb + r
    return best, cx, cy, r
end

# --- TM style: short radial ticks (bars) ---
function add_tm_ticks!(
    plt, x::AbstractVector, y::AbstractVector, rotation::AbstractVector, r::Real;
    indices = eachindex(x),
    every::Int = 1,
    len_frac::Real = 0.35,     # length of the bar relative to r
    inset_frac::Real = 0.10,   # how far inside the rim the bar starts
    lw::Union{Real,Symbol} = :auto,
    colors = (:blue, :red, :yellow, :green),
    alpha::Real = 0.95,
)
    # adaptive linewidth (in points)
    L = length(indices)
    lw_pts = lw === :auto ? max(0.8, 2.8 - 1.2*log10(max(L, 10))) : Float64(lw)

    @inbounds for k in 1:every:length(indices)
        i = indices[k]
        xi = x[i]; yi = y[i]
        rot0 = rotation[i]

        for (j, adddeg) in enumerate((0.0, 90.0, 180.0, 270.0))
            θ = (rot0 + adddeg) * π/180
            c, s = cos(θ), sin(θ)

            r0 = r*(1 - inset_frac)
            r1 = r*(1 - inset_frac + len_frac)
            x0 = xi + r0*c; y0 = yi + r0*s
            x1 = xi + r1*c; y1 = yi + r1*s

            plot!(plt, [x0, x1], [y0, y1];
                  seriestype=:path, color=colors[j], alpha=alpha,
                  linewidth=lw_pts, linecap=:round, label=false)
        end
    end
    return plt
end

# --- TM style: tiny triangles (paddles) ---  (VECTOR-ONLY; no Shape)
function add_tm_triangles!(
    plt, x::AbstractVector, y::AbstractVector, rotation::AbstractVector, r::Real;
    indices = eachindex(x),
    every::Int = 1,
    size_frac::Real = 0.28,  # overall marker size control
    colors = (:blue, :red, :yellow, :green),
    alpha::Real = 0.95,
)
    # geometry factors (relative to r)
    tip_off   = 1.05
    base_rad  = 0.85
    base_half = 0.16
    # mild scaling around default 0.28
    tip_off   *= (1 + 0.6*(size_frac-0.28)/0.28)
    base_rad  *= (1 - 0.4*(size_frac-0.28)/0.28)
    base_half *= (1 + 0.8*(size_frac-0.28)/0.28)

    @inbounds for k in 1:every:length(indices)
        i = indices[k]
        xi = x[i]; yi = y[i]
        rot0 = rotation[i]

        # skip any bad entries
        (isfinite(xi) && isfinite(yi) && isfinite(rot0)) || continue

        for (j, adddeg) in enumerate((0.0, 90.0, 180.0, 270.0))
            θ = (rot0 + adddeg) * π/180
            c, s = cos(θ), sin(θ)
            ux, uy =  c,  s
            tx, ty = -s,  c

            tipx  = xi + (r*tip_off)*ux
            tipy  = yi + (r*tip_off)*uy
            basex = xi + (r*base_rad)*ux
            basey = yi + (r*base_rad)*uy
            v1x   = basex + (r*base_half)*tx
            v1y   = basey + (r*base_half)*ty
            v2x   = basex - (r*base_half)*tx
            v2y   = basey - (r*base_half)*ty

            # build as VECTORS (not tuples) — avoids Float64 series error
            xs = Float64[tipx, v1x, v2x]
            ys = Float64[tipy, v1y, v2y]

            # final guard (just in case)
            if all(isfinite, xs) && all(isfinite, ys)
                plot!(plt, xs, ys;
                      seriestype = :shape,
                      fillcolor  = colors[j],
                      fillalpha  = alpha,
                      linecolor  = :transparent,
                      label      = false)
            end
        end
    end
    return plt
end

# Batched backbone lines (your function; kept here for locality)
function _draw_edges_batched!(plt, x::AbstractVector, y::AbstractVector,
                              edges::Vector{Tuple{Int,Int}};
                              lw::Real=1.5, alpha::Real=0.9, color=:black)
    n = length(edges)
    xs = Vector{Float64}(undef, 3n)
    ys = Vector{Float64}(undef, 3n)
    k = 1
    @inbounds for (u,v) in edges
        xs[k] = x[u]; ys[k] = y[u]; k += 1
        xs[k] = x[v]; ys[k] = y[v]; k += 1
        xs[k] = NaN;  ys[k] = NaN;  k += 1
    end
    plot!(plt, xs, ys; seriestype=:path, lw=lw, alpha=alpha, color=color, label=false)
    return plt
end

"""
    draw_backbone_with_roles!(plt, x, y, edges;
        turn_thresh_deg=30,
        role_ms=6,
        contour_color=:black, contour_width=1.5)

Draw backbone lines + role overlays using same contour/marker style you use elsewhere.
"""
function draw_backbone_with_roles!(plt, x::AbstractVector, y::AbstractVector,
                                   edges::Vector{Tuple{Int,Int}};
                                   turn_thresh_deg::Real=30)

    N = length(x)

    # --- base size scaling ---
    base = clamp(8000.0 / N, 0.4, 6.0)   # 2000 → ~1, 8000 → ~0.5, floor at 0.4, cap at 6

    ms_turn = base
    ms_leaf = base * 1.5
    ms_junc = base * 1.35

    # 1) Backbone lines
    _draw_edges_batched!(plt, x, y, edges; lw=2, alpha=0.9, color=:black)

    # 2) Roles
    roles, _, _ = classify_backbone_nodes(x, y, edges; turn_thresh_deg=turn_thresh_deg)
    leaf_idx  = findall(r -> r === :leaf, roles)
    junc_idx  = findall(r -> r === :junction, roles)
    turn_idx  = findall(r -> r === :turning, roles)

    # 3) Overlays — fill only, no contour
    if !isempty(leaf_idx)
        scatter!(plt, x[leaf_idx], y[leaf_idx];
            marker=:circle, ms=ms_leaf,
            markercolor=:cyan, markerstrokewidth=0, label=false)
    end
    if !isempty(junc_idx)
        scatter!(plt, x[junc_idx], y[junc_idx];
            marker=:circle, ms=ms_junc,
            markercolor=:forestgreen, markerstrokewidth=0, label=false)
    end
    if !isempty(turn_idx)
        scatter!(plt, x[turn_idx], y[turn_idx];
            marker=:square, ms=ms_turn,
            markercolor=:orange, markerstrokewidth=0, label=false)
    end

    return plt
end



# Map vertex -> degree
degree_map(edges::Vector{Tuple{Int,Int}}) = begin
    deg = Dict{Int,Int}()
    @inbounds for (u,v) in edges
        deg[u] = get(deg, u, 0) + 1
        deg[v] = get(deg, v, 0) + 1
    end
    deg
end

# Angle deviation at deg==2 nodes (deviation from straight, in degrees)
function _turn_deviation_deg(i::Int, x, y, nbrs::Vector{Vector{Int}})
    n1, n2 = nbrs[i][1], nbrs[i][2]
    v1x, v1y = x[n1] - x[i], y[n1] - y[i]
    v2x, v2y = x[n2] - x[i] - 0.0, y[n2] - y[i] - 0.0
    l1 = hypot(v1x, v1y); l2 = hypot(v2x, v2y)
    if l1 == 0 || l2 == 0
        return 180.0  # treat as a hard bend for safety
    end
    dot12 = clamp((v1x*v2x + v1y*v2y) / (l1*l2), -1.0, 1.0)
    θ = acosd(dot12)                 # 0..180 (interior)
    abs(180.0 - θ)                   # deviation from straight
end

# Build adjacency
function _adjacency(N::Int, edges::Vector{Tuple{Int,Int}})
    nbrs = [Int[] for _ in 1:N]
    @inbounds for (u,v) in edges
        push!(nbrs[u], v)
        push!(nbrs[v], u)
    end
    nbrs
end

# Classify nodes
function classify_backbone_nodes(x::AbstractVector, y::AbstractVector,
                                 edges::Vector{Tuple{Int,Int}};
                                 turn_thresh_deg::Real = 15)
    N = length(x)
    nbrs = _adjacency(N, edges)
    deg = [length(nbrs[i]) for i in 1:N]
    roles = fill(:isolated, N)

    @inbounds for i in 1:N
        if deg[i] == 0
            roles[i] = :isolated
        elseif deg[i] == 1
            roles[i] = :leaf
        elseif deg[i] >= 3
            roles[i] = :junction
        else
            # deg == 2 → turning vs interim
            dev = _turn_deviation_deg(i, x, y, nbrs)
            roles[i] = (dev >= turn_thresh_deg) ? :turning : :interim_straight
        end
    end

    return roles, deg, nbrs
end


# -----------------------------
# Utilities
# -----------------------------
# replace your density_curve with this:

function density_curve(data; label="", xlabel="", npts=512)
    # accept (xs, ys) or (idxs, vals) tuples by taking the 2nd thing,
    # otherwise assume it's a vector of numbers
    if data isa Tuple && length(data) == 2
        data = data[2]
    end
    d = collect(skipmissing(data))           # materialize: no SkipMissing iter
    isempty(d) && return plot(title="(no data)"; legend=false)

    kd = kde(d)                              # avoid npoints kw to dodge MethodError
    plot(kd.x, kd.density; seriestype=:path, lw=2,
         label=label, xlabel=xlabel, ylabel="density", legend=:topright)
end


# 2) Segment length
segment_length(x, y, seg::Vector{Int}) = sum(hypot(x[seg[i+1]]-x[seg[i]], y[seg[i+1]]-y[seg[i]]) for i in 1:length(seg)-1)

# 3) Angle deviation from straight (180°) at a triple (a,b,c)
@inline function turn_dev_deg(ax,ay,bx,by,cx,cy)
    v1x, v1y = ax-bx, ay-by
    v2x, v2y = cx-bx, cy-by
    l1 = hypot(v1x,v1y); l2 = hypot(v2x,v2y)
    (l1==0 || l2==0) && return 180.0
    cosθ = clamp((v1x*v2x + v1y*v2y)/(l1*l2), -1.0, 1.0)
    θ = acosd(cosθ)              # 0..180
    abs(180.0 - θ)               # deviation from straight
end

# 4) Per-segment curvature stats: (total_turn, mean_turn, curvature_density)
function segment_curvature_stats(x, y, seg::Vector{Int})
    n = length(seg)
    if n < 3
        return (0.0, 0.0, 0.0)
    end
    devs = Float64[]
    @inbounds for i in 2:n-1
        a, b, c = seg[i-1], seg[i], seg[i+1]
        push!(devs, turn_dev_deg(x[a],y[a], x[b],y[b], x[c],y[c]))
    end
    total = sum(devs) * (π/180)  # convert deg -> rad for totals
    meanv = mean(devs) * (π/180)
    len   = segment_length(x,y,seg)
    dens  = (len > 0 ? total/len : 0.0)  # rad / unit length
    return (total, meanv, dens)
end

# -----------------------------
# Graph-wide features
# -----------------------------
# B) Turning-angle deviations across the whole graph (degree-2 interior turns)
function graph_turn_deviation_deg(x, y, edges::Vector{Tuple{Int,Int}})
    N = length(x)
    nbrs = [Int[] for _ in 1:N]
    @inbounds for (u,v) in edges
        push!(nbrs[u], v); push!(nbrs[v], u)
    end
    devs = Float64[]
    @inbounds for b in 1:N
        if length(nbrs[b]) == 2
            a, c = nbrs[b][1], nbrs[b][2]
            push!(devs, turn_dev_deg(x[a],y[a], x[b],y[b], x[c],y[c]))
        end
    end
    return devs
end

# C) Per-segment curvature density vector
function per_segment_curvature_density(x,y, segments)
    curv_stats = [segment_curvature_stats(x,y,seg) for seg in segments]
    getindex.(curv_stats, 3)  # dens
end

"""
    segment_lengths(x, y, segments) -> Vector{Float64}

Compute the Euclidean length of each segment.
- `x`, `y`: node coordinates
- `segments`: list of node-index sequences (each segment is a polyline)

Returns a vector of lengths (one per segment).
"""
function segment_lengths(x::AbstractVector, y::AbstractVector,
                         segments::Vector{Vector{Int}})
    lengths = Float64[]
    for seg in segments
        L = 0.0
        for i in 1:length(seg)-1
            dx = x[seg[i+1]] - x[seg[i]]
            dy = y[seg[i+1]] - y[seg[i]]
            L += hypot(dx, dy)
        end
        push!(lengths, L)
    end
    return lengths
end



# -----------------------------
# Plotters (smooth curves)
# -----------------------------
plot_segment_length_density(x, y, segments) = begin
    lens = segment_lengths_safe(x, y, segments)
    plt  = density_curve(lens; label="segment length", xlabel="length (plot units)")
    return lens, plt
end

# keep this (vector-of-edges style):
plot_turning_angle_density(x, y, edges) =
    density_curve(graph_turn_deviation_deg(x, y, edges);
                  label="turning angle dev.", xlabel="degrees")

# add this overload (vector-of-devs style):
plot_turning_angle_density(devs::AbstractVector) =
    density_curve(devs; label="turning angle dev.", xlabel="degrees")

plot_curvature_density(x,y, segments) =
    density_curve(per_segment_curvature_density(x,y,segments); label="curvature density", xlabel="rad / unit length")

function normalize_edges!(edges::Vector{Tuple{Int,Int}})
    @inbounds for k in eachindex(edges)
        u,v = edges[k]
        u == v && error("self-loop at $(u)")
        edges[k] = (min(u,v), max(u,v))
    end
    unique!(edges)
    return edges
end

# keep at top-level near your other helpers
function _dedup_run(seg::Vector{Int})
    out = Int[]
    last = typemin(Int)
    @inbounds for i in seg
        if i != last
            push!(out, i); last = i
        end
    end
    out
end

function segment_lengths_safe(x::AbstractVector, y::AbstractVector,
                              segments::Vector{Vector{Int}}; min_len=1e-9)
    lens = Float64[]
    @inbounds for seg in segments
        s = _dedup_run(seg)
        if length(s) ≥ 2
            L = 0.0
            for k in 1:length(s)-1
                i, j = s[k], s[k+1]
                L += hypot(x[j]-x[i], y[j]-y[i])
            end
            (isfinite(L) && L > min_len) && push!(lens, L)
        end
    end
    return lens
end


# --- convenience: run face-finder on arbitrary x,y without touching state ---
# Requires the Faces module we wrote earlier to be in scope.
struct XYEdges
    x_coords::Vector{Float64}
    y_coords::Vector{Float64}
    edges::Vector{Tuple{Int,Int}}
end

compute_enclosed_faces_xy(x::AbstractVector{<:Real},
                          y::AbstractVector{<:Real},
                          edges::Vector{Tuple{Int,Int}};
                          area_floor::Float64=0.0,
                          drop_outer::Bool=true,
                          normalize_orientation::Bool=true,
                          return_abs::Bool=true) =
    Faces.compute_enclosed_faces(XYEdges(collect(x), collect(y), edges);
        area_floor=area_floor,
        drop_outer=drop_outer,
        normalize_orientation=normalize_orientation,
        return_abs=return_abs)

# --- plot faces onto an existing Plots.jl plot ---
function draw_faces_on!(plt, x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                        edges::Vector{Tuple{Int,Int}};
                        inputstate = nothing,
                        min_area::Float64 = 5.0,    # filter tiny faces by |area|
                        max_faces::Int = 5_000,     # cap draw count for speed
                        face_alpha::Float64 = 0.35,
                        outline_color = :gray20,
                        outline_lw::Real = 0.5,
                        color_by::Symbol = :area,   # :area or :uniform
                        show_labels::Bool = false,
                        label_fontsize::Int = 7)


    println(typeof(inputstate))
    res = Faces.compute_enclosed_faces(inputstate; area_floor=min_area, drop_outer=true)

    faces = res.faces
    areas = res.areas                   # already abs if return_abs=true
    println(length(areas))
    # keep faces by area and limit count
    keep = findall(a -> a ≥ min_area, areas)
    isempty(keep) && return plt
    ord = sort(keep; by=i->areas[i], rev=true)
    ord = ord[1:min(length(ord), max_faces)]

    # min-max normalize areas for coloring
    amin, amax = extrema(areas[ord]); rng = max(amax-amin, eps())

    # centroid helper for labels
    function poly_centroid(vs::Vector{Int})
        n = length(vs)
        acc = 0.0; cx = 0.0; cy = 0.0
        @inbounds for i in 1:n
            j = (i == n) ? 1 : i+1
            xi, yi = x[vs[i]], y[vs[i]]
            xj, yj = x[vs[j]], y[vs[j]]
            cr = xi*yj - xj*yi
            acc += cr
            cx += (xi + xj) * cr
            cy += (yi + yj) * cr
        end
        if abs(acc) < 1e-15
            return (mean(x[vs]), mean(y[vs]))
        else
            acc2 = acc * 3.0
            return (cx/acc2, cy/acc2)
        end
    end

    for i in ord
        vs = faces[i]
        xs = x[vs]; ys = y[vs]

        # color
        c = :dodgerblue
        if color_by === :area
            t = (areas[i]-amin)/rng
            # simple two-stop blue gradient
            r = (1-t)*0.8 + t*0.0
            g = (1-t)*0.9 + t*0.2
            b = (1-t)*1.0 + t*0.8
            c = RGB(1,0,0)
        elseif color_by === :uniform
            c = RGBA(0.2,0.6,0.9,1.0)
        end

        # filled polygon (close it)
        plot!(plt, vcat(xs, xs[1]), vcat(ys, ys[1]);
              seriestype=:shape, fillalpha=face_alpha, fillcolor=c,
              linecolor=outline_color, lw=outline_lw)

        if show_labels
            cx, cy = poly_centroid(vs)
            annotate!(plt, (cx, cy, text(string(round(areas[i]; digits=2)), label_fontsize, :black)))
        end
    end
    return plt
end


# -------------------------
# Main LOD plotter
# -------------------------
"""
plot_monomers_lod(x, y; rotation, boxSize, monomer_radius, ...)

Features:
- Green-filled circles, red contour
- Spin shown as a green rim segment pointing to `rotation[i]` (deg)
- Auto units nm/μm/mm using NM_PER_DATA map
- LOD: detail/medium/massive
"""
function plot_monomers_lod(
    x::AbstractVector, y::AbstractVector;
    backbones = nothing,
    inputstate = nothing,
    rotation::Union{AbstractVector,Nothing}=nothing,
    boxSize::Union{Real,Nothing}=nothing,
    monomer_radius::Union{Nothing,Real}=nothing,
    mode=:min, scale=:identity, target_max=nothing,
    tick_step::Union{Real,Symbol,Nothing}=:auto,
    xlim=nothing, ylim=nothing,
    show_grid::Bool=false, show_orientation::Bool=true,
    lod::Symbol=:auto, sample_cap_points::Int=100_000,
    bins::Tuple{Int,Int}=(600,600),
    orient_every::Union{Int,Symbol}=:auto,
    draw_contacts::Union{Bool,Symbol}=:auto,
    contact_scale::Real=1.02,
    contact_max_lines::Int=200_000,
    orient_len::Union{Real,Symbol}=:auto,
    marker_size::Union{Real,Symbol}=:auto,
    figsize::Tuple{Int,Int}=(1200,1200),
    dpi::Int=200,
    # physical units
    real_scale_nm::Real = NM_PER_DATA,
    unit::Union{Symbol,String} = :auto,
    save_path=nothing,
    # NEW: outline + TM style controls
    show_contour::Bool=true,
    tm_style::Union{Symbol,String}=:auto,
    tm_len_frac::Real=0.35,
    tm_inset_frac::Real=0.10,
    tm_lw::Union{Real,Symbol}=:auto,
    tm_size_frac::Real=0.28,
    # in the function signature, add:
    edges::Union{Nothing,Vector{Tuple{Int,Int}}}=nothing,

)
    N = length(x)
    if isempty(x) || isempty(y)
        @warn "No points to plot"
        return
    end

    # 1) transform once
    xT, yT, s = transform_coords(x, y; mode=mode, scale=scale, target_max=target_max)
    normalize_edges!(edges)  # NEW

    # 2) physical units → chosen display unit
    x_nm = xT .* real_scale_nm
    y_nm = yT .* real_scale_nm
    r_nm = monomer_radius === nothing ? nothing : Float64(monomer_radius) * s * real_scale_nm

    spanx_nm = maximum(x_nm) - minimum(x_nm)
    spany_nm = maximum(y_nm) - minimum(y_nm)
    span_nm  = max(spanx_nm, spany_nm)
    u   = unit === :auto ? pick_unit(span_nm) : Symbol(unit)
    div = unit_divisor(u)

    xT = x_nm ./ div
    yT = y_nm ./ div
    r_plot_unit = r_nm === nothing ? nothing : r_nm / div

    # 3) LOD
    if lod == :auto
        lod = N ≤ 100 ? :detail : (N ≤ 800 ? :medium : :massive)
    end

    # 4) limits + ticks
    spanx = maximum(xT) - minimum(xT)
    spany = maximum(yT) - minimum(yT)
    pad = r_plot_unit !== nothing ? 1.15 * r_plot_unit : 0.02 * max(spanx, spany)

    if xlim === nothing || ylim === nothing
        xlim = (minimum(xT) - pad, maximum(xT) + pad)
        ylim = (minimum(yT) - pad, maximum(yT) + pad)
    else
        xlim = (first(xlim) - pad, last(xlim) + pad)
        ylim = (first(ylim) - pad, last(ylim) + pad)
    end

    xtick_vals = nothing
    ytick_vals = nothing
    if tick_step === :auto
        step = nice_tick_step(max(last(xlim)-first(xlim), last(ylim)-first(ylim)))
        x0 = first(xlim) - mod(first(xlim), step)
        y0 = first(ylim) - mod(first(ylim), step)
        xtick_vals = collect(x0:step:last(xlim))
        ytick_vals = collect(y0:step:last(ylim))
    elseif tick_step !== nothing
        step = Float64(tick_step)
        xtick_vals = collect(first(xlim):step:last(xlim))
        ytick_vals = collect(first(ylim):step:last(ylim))
    end

    # 5) marker size
    ms_from_radius = nothing
    if r_plot_unit !== nothing && lod != :massive
        data_diam_points = function(r_data, xlim, ylim, figsize, dpi)
            wpx, hpx = figsize
            spanx_ = last(xlim) - first(xlim)
            spany_ = last(ylim) - first(ylim)
            spanu  = max(spanx_, spany_)
            px_per_unit = min(wpx, hpx) / spanu
            diam_px = 2 * r_data * px_per_unit
            return diam_px * 72 * 1.35 / dpi
        end
        ms_from_radius = data_diam_points(r_plot_unit, xlim, ylim, figsize, dpi) * 1.05
    end

    px_scale = max(figsize...) / 800
    ms_auto  = (lod == :detail ? 8.0 : (lod == :medium ? 4.0 : 0.0)) * px_scale
    ms_final = marker_size === :auto ? (ms_from_radius === nothing ? ms_auto : ms_from_radius) : Float64(marker_size)

    # 6) orientation sampling
    target_arrows = lod == :detail ? 200 : 120
    oe_auto = max(1, Int(ceil(N / target_arrows)))
    oe = orient_every === :auto ? oe_auto : Int(orient_every)

    # 7) plot
    plt = nothing
    plt2 = nothing 
    if lod == :massive        # Base figure with axes/ticks, no monomer markers
        plt = plot(
            aspect_ratio=:equal, legend=false, dpi=dpi, size=figsize,
            xlim=xlim, ylim=ylim, xticks=xtick_vals, yticks=ytick_vals,
            grid=false
        )
        plt2 = plot(
            aspect_ratio=:equal, legend=false, dpi=dpi, size=figsize,
            xlim=xlim, ylim=ylim, xticks=xtick_vals, yticks=ytick_vals,
            grid=false
        )
        plt_backbone = plot(
            aspect_ratio=:equal, legend=false, dpi=dpi, size=figsize,
            xlim=xlim, ylim=ylim, xticks=xtick_vals, yticks=ytick_vals,
            grid=false
        )

        if edges !== nothing && !isempty(edges)
            
            # draw backbone ONLY (no points) using transformed coords (xT,yT)
            _draw_edges_batched!(plt, xT, yT, edges; lw=0.7, alpha=0.85, color=:black)
            # scatter!([x[1] for x in strand_starts], [x[2] for x in strand_starts]; m=:xcross, ms=6, mc=:red, label="strand starts")

        else
            @warn "plot_monomers_lod: edges not provided; massive mode wants precomputed backbone"
        end

        # overlay faces on the SAME transformed coords
        draw_faces_on!(plt, xT, yT, edges;
            inputstate = inputstate,
            min_area = 1.0,       # tune to your scale to hide slivers
            max_faces = 3000,      # keep plots responsive
            face_alpha =1.0,
            color_by = :area,
            show_labels = false
        )


        for cid in sort(collect(keys(backbones)))
            path = backbones[cid]
            isempty(path) && continue

            if length(path) >= 2
                # polyline for this island
                plot!(plt2, xT[path], yT[path], lw=0.7, alpha=0.95, label = "island $cid")
                # endpoints
                scatter!(plt2, [xT[first(path)]], [yT[first(path)]], ms=.5, label="")
                scatter!(plt2, [xT[last(path)]],  [yT[last(path)]],  ms=.5, label="")
            else
                # singleton island (no edges) — draw a dot
                # scatter!(plt, [xT[path[1]]], [yT[path[1]]], ms=.5, label = "island $cid (singleton)")
            end
        end

        # skip the rest of the body for massive mode
        xlabel!(plt, "Distance ($(unit_label(u)))")
        ylabel!(plt, "Distance ($(unit_label(u)))")
        if save_path !== nothing
            save_path_strand   = replace(save_path, "placement" => "strand")
            save_path_backbone = replace(save_path, "placement" => "backbone")

            savefig(plt, save_path_strand)
            savefig(plt2, save_path_backbone)
            println("Plot saved to: $save_path")
        end
        



        draw_backbone_with_roles!(plt_backbone, xT, yT, edges;
            turn_thresh_deg = 30.0,
        )
        display(plt)
        display(plt_backbone)
        out_dir = joinpath(pwd(), "plots", "tmp")

        if save_path !== nothing
            savefig(plt_backbone, out_dir*"\\bb.png")
            println("Plot saved to: $save_path")
        end

        # segments, endpoints, junctions = segments_from_backbone(xT, yT, edges)
        # segments, endpoints, junctions = segments_from_backbone_cc(xT, yT, edges)

        # seg_lengths, plt_len = plot_segment_length_density(xT, yT, segments)
        # plt_turn = plot_turning_angle_density(xT, yT, edges)
        # plt_curv = plot_curvature_density(xT, yT, segments)
        # # Compute in nanometers
        # lens_nm = segment_lengths_safe(x_nm, y_nm, segments)


        # # KDE bandwidth tuned for nm
        # kd = kde(lens_nm; bandwidth=1.5)
        # mask = kd.x .>= 0.0
        
        # plt_len = plot(kd.x[mask], kd.density[mask];
        #     xlabel="length (nm)", ylabel="density",
        #     title="Segment length distribution (KDE)",
        #     lw=2, color=:blue, legend=false)
        # histogram!(lens_nm, normalize=true, alpha=0.3, bins=50, color=:blue)


        # display(plt_len)
        # println("max segment length: $(maximum(lens_nm)) nm")

        # # optional save
        # savefig(plt_len, replace(save_path, ".png" => "_seglen_kde.png"))

        # savefig(plt_turn, replace(save_path, ".png" => "_turn_kde.png"))
        # savefig(plt_curv, replace(save_path, ".png" => "_curvdens_kde.png"))



        return plt

    else
        idx = 1:length(xT)
        if lod == :medium && N > sample_cap_points
            stride = ceil(Int, N / sample_cap_points)
            idx = 1:stride:length(xT)
        end

        # contour toggle
        contour_color = show_contour ? :black : :transparent
        contour_width = show_contour ? 2.0 : 0.0

        plt = scatter(xT[idx], yT[idx];
            marker=:circle, ms=ms_final,
            markercolor=:white,
            markerstrokecolor=contour_color,
            markerstrokewidth=contour_width,
            linecolor=:transparent, aspect_ratio=:equal, legend=false, grid=false,
            dpi=dpi, size=figsize,
            xlim=xlim, ylim=ylim, xticks=xtick_vals, yticks=ytick_vals)

        # --- Backbone (finder + raw draw) -------------------------
        do_contacts = (draw_contacts === :auto ? (lod != :massive) : draw_contacts) && (r_plot_unit !== nothing)
        edges = Tuple{Int,Int}[]
        if do_contacts
            # Use FULL xT/yT for geometry; don't subsample with idx
            edges = find_contact_edges(xT, yT, r_plot_unit;
                                    xlim=xlim, ylim=ylim,
                                    contact_scale=contact_scale,
                                    max_edges=contact_max_lines)
            draw_backbone!(plt, xT, yT, edges; lc=:black, lw=1.5, alpha=0.9)  # f2
        end



        # --- Standalone backbone figure with colored segments + vertices -------------
        plt_backbone = plot(; size=(1600,1600), aspect_ratio=:equal, legend=false, xlim=xlim, ylim=ylim)
        if r_plot_unit !== nothing
            edges_back = isempty(edges) ? find_contact_edges(xT, yT, r_plot_unit;
                                                            xlim=xlim, ylim=ylim,
                                                            contact_scale=contact_scale,
                                                            max_edges=contact_max_lines) : edges
            # base backbone
            draw_backbone!(plt_backbone, xT, yT, edges_back; lc=:black, lw=2, alpha=0.9)

            # DEGREE-BASED FEATURES (f1)
            # feat = backbone_features(xT, yT, edges_back)
            # compare_edges(state, edges)

            # 💡 CURVATURE REFINEMENT (split at turning points)
            # 1) detect + split with looser params
            # segments2, tverts = refine_segments_by_curvature(
            #     xT, yT, feat.segments;
            #     angle_thresh_deg = 30.0,   # ↓ from 60
            #     min_gap = 1                # ↓ from 2
            # )

            

            tverts = turning_vertices_from_edges(xT, yT, edges_back; angle_thresh_deg=30.0)
            @info "Graph-turns" n_turns = length(tverts)

            # A) Segment length: smooth density
            seg_lengths, plt_len = plot_segment_length_density(
                xT, yT, feat.segments;
            )

            # D) Turning angles: smooth density (graph-based)
            _, dev_deg_all = graph_turn_deviation_deg(xT, yT, edges_back)
            plt_ta = plot_turning_angle_density(dev_deg_all;
                save_path = (save_path === nothing ? nothing : replace(save_path, ".png" => "_turning_angles_density.png"))
            )

            ###############
            # (B) Segment curvature stats (per segment)
            ###############
            curv_stats = [ segment_curvature_stats(xT, yT, seg) for seg in feat.segments ]
            totals  = getindex.(curv_stats, 1)
            means   = getindex.(curv_stats, 2)
            dens    = getindex.(curv_stats, 3)

            plt_curv_total = histogram(totals; bins=30, xlabel="total |turn| (rad)", ylabel="freq",
                                    title="Per-segment total curvature", legend=false)
            plt_curv_mean  = histogram(means;  bins=30, xlabel="mean |turn| per interior (rad)", ylabel="freq",
                                    title="Per-segment mean curvature", legend=false)
            plt_curv_dens  = histogram(dens;   bins=30, xlabel="curvature density (rad / unit length)", ylabel="freq",
                                    title="Per-segment curvature density", legend=false)

            if save_path !== nothing
                savefig(plt_curv_total, replace(save_path, ".png" => "_curv_total_hist.png"))
                savefig(plt_curv_mean,  replace(save_path, ".png" => "_curv_mean_hist.png"))
                savefig(plt_curv_dens,  replace(save_path, ".png" => "_curv_density_hist.png"))
            end

            ###############
            # (C) Junction and leaf vertices
            ###############
            deg = degree_map(edges_back)
            leaf_vertices    = [v for (v,d) in deg if d == 1]
            junction_vertices = [v for (v,d) in deg if d ≥ 3]
            @info "Graph degrees" n_leaf=length(leaf_vertices) n_junction=length(junction_vertices)

            # Optional: branching (degree) frequency at junctions
            branching_counts = [deg[v] for v in junction_vertices]
            if !isempty(branching_counts)
                m = maximum(branching_counts)
                plt_branches = histogram(branching_counts; bins=1:1:m+1,
                    xticks=(1:m, string.(1:m)),
                    xlabel="degree at junction", ylabel="freq", title="Junction branching degree", legend=false)
                if save_path !== nothing
                    savefig(plt_branches, replace(save_path, ".png" => "_junction_branching_hist.png"))
                end
            end


            # 2) overlay: PASS turning_vertices directly (so they’re always shown)
            overlay_segments_vertices!(
                plt_backbone, xT, yT,
                (; segments = feat.segments,
                endpoints = feat.endpoints,
                junctions = feat.junctions);
                segment_mode = :cmap,
                v_ms = 8,
                turning_vertices = tverts,       # <— NEW: wire them in
                turning_color = :orange,
                turning_ms = 4,
                turning_marker = :square
            )
            

            if save_path !== nothing
                savefig(plt_backbone, replace(save_path, ".png" => "_backbone.png"))
            end
        end

        # compass legend
        begin
            _, cx, cy, rcomp = _choose_compass_corner(xT, yT, xlim, ylim;
                                                      size_frac=0.12, margin_frac=0.05,
                                                      sample_cap=sample_cap_points)
            θ = range(0, 2π; length=200)
            plot!(plt, cx .+ rcomp .* cos.(θ), cy .+ rcomp .* sin.(θ);
                seriestype=:path, color=:black, linewidth=1.5,
                fillcolor=RGBA(0,0,0,0), label=false)

            dirs = [
                (0.0,        :blue,   "0/360"),
                (π/2,        :red,    "90"),
                (π,          :yellow, "180"),
                (3π/2,       :green,  "270"),
            ]
            for (ang, col, lab) in dirs
                x1, y1 = cx + 0.65*rcomp*cos(ang), cy + 0.65*rcomp*sin(ang)
                x2, y2 = cx + 1.00*rcomp*cos(ang), cy + 1.00*rcomp*sin(ang)
                plot!(plt, [x1, x2], [y1, y2]; color=col, linewidth=2.0, linecap=:round, label=false)
                xt, yt = cx + 1.22*rcomp*cos(ang), cy + 1.22*rcomp*sin(ang)
                annotate!(plt, xt, yt, text(lab, 9, col, :center))
            end
        end

        # TM markers (skip if tm_style says "none")
        tstyle = Symbol(tm_style)
        if show_orientation && rotation !== nothing && !isempty(rotation) &&
        !(tstyle in (:nothing, :none))
            maxidx = maximum(collect(idx))
            if length(rotation) ≥ maxidx
                vis_ids = collect(idx)
                if r_plot_unit !== nothing
                    style = tstyle == :auto ? (show_contour ? :arcs : :ticks) : tstyle
                    if style == :arcs
                        add_orientation_arcs!(plt, xT, yT, rotation, r_plot_unit;
                                              indices=vis_ids, every=oe,
                                              arc_span_deg=25.0, thickness_frac=0.45,
                                              color=:blue,  alpha=0.95)
                        add_orientation_arcs!(plt, xT, yT, rotation .+ 90.0, r_plot_unit;
                                              indices=vis_ids, every=oe,
                                              arc_span_deg=25.0, thickness_frac=0.45,
                                              color=:red,   alpha=0.95)
                        add_orientation_arcs!(plt, xT, yT, rotation .+ 180.0, r_plot_unit;
                                              indices=vis_ids, every=oe,
                                              arc_span_deg=25.0, thickness_frac=0.45,
                                              color=:yellow,alpha=0.95)
                        add_orientation_arcs!(plt, xT, yT, rotation .+ 270.0, r_plot_unit;
                                              indices=vis_ids, every=oe,
                                              arc_span_deg=25.0, thickness_frac=0.45,
                                              color=:green, alpha=0.95)
                    elseif style == :ticks
                        add_tm_ticks!(plt, xT, yT, rotation, r_plot_unit;
                                      indices=vis_ids, every=oe,
                                      len_frac=tm_len_frac, inset_frac=tm_inset_frac,
                                      lw=tm_lw, colors=(:blue,:red,:yellow,:green), alpha=0.95)
                    elseif style == :triangles
                        add_tm_triangles!(plt, xT, yT, rotation, r_plot_unit;
                                          indices=vis_ids, every=oe,
                                          size_frac=tm_size_frac,
                                          colors=(:blue,:red,:yellow,:green), alpha=0.95)
                    else
                        @warn "tm_style=$style not recognized; falling back to :ticks"
                        add_tm_ticks!(plt, xT, yT, rotation, r_plot_unit;
                                      indices=vis_ids, every=oe,
                                      len_frac=tm_len_frac, inset_frac=tm_inset_frac,
                                      lw=tm_lw, colors=(:blue,:red,:yellow,:green), alpha=0.95)
                    end
                else
                    oidx = vis_ids[1:oe:length(vis_ids)]
                    spanu = max(spanx, spany)
                    alen = 0.03 * spanu
                    u = cosd.(rotation[oidx]); v = sind.(rotation[oidx])
                    quiver!(plt, xT[oidx], yT[oidx], quiver=(u .* alen, v .* alen);
                            lw=0.5, alpha=0.8, linecolor=:black, label=false)
                end
            end
        end
    end

    # optional grid overlay
    if show_grid && boxSize !== nothing
        scaled_box = Float64(boxSize) * s * (real_scale_nm / div)
        x_end, y_end = last(xlim), last(ylim)
        for xg in 0:scaled_box:x_end
            plot!(plt, [xg, xg], [first(ylim), y_end], lw=0.5, alpha=0.3, linecolor=:gray, label=false)
        end
        for yg in 0:scaled_box:y_end
            plot!(plt, [first(xlim), x_end], [yg, yg], lw=0.5, alpha=0.3, linecolor=:gray, label=false)
        end
    end

    xlabel!(plt, "Distance ($(unit_label(u)))")
    ylabel!(plt, "Distance ($(unit_label(u)))")

    if save_path !== nothing
        savefig(plt, save_path)
        println("Plot saved to: $save_path")
    end

    return plt
end


function plot_W_vs_K_max(W::AbstractMatrix, K::AbstractMatrix;
                         add_angle_ticks::Bool=true,
                         save_path::Union{Nothing,String}=nothing)

    @assert size(W) == size(K) "W and K must have same shape (e.g., 72×72)."
    ZED = cgrad(:roma, 10, categorical=true, scale=:exp)

    # --- Normalize both W and K by their own max ---
    Wn = W ./ max(maximum(W), eps())
    Kn = K ./ max(maximum(K), eps())

    clims = (0.0, 1.0)  # shared colorbar limits

    pltW = heatmap(Wn;
        aspect_ratio=:equal,
        color=ZED,
        clims=clims,
        colorbar_title="normalized (0–1)",
        title="Input connection probability (W)")

    pltK = heatmap(Kn;
        aspect_ratio=:equal,
        color=ZED,
        clims=clims,
        colorbar_title="normalized (0–1)",
        title="Simulated connections (K)")

    # --- Degree ticks (0°,90°,180°,270°) for 72 bins (5° each) ---
    if add_angle_ticks && size(W,1) == 72 && size(W,2) == 72
        ticks = 1:18:72
        labs  = ["0°","90°","180°","270°"]
        for p in (pltW, pltK)
            xticks!(p, (ticks, labs))
            yticks!(p, (ticks, labs))
        end
    end

    fig = plot(pltW, pltK; layout=(1,2), size=(1200,520), dpi=200)
    if save_path !== nothing
        savefig(fig, save_path)
        println("Saved → $save_path")
    end
    return fig
end


# -------------------------
# Public entry that uses state/config
# -------------------------

function generate_plots(state::AbstractState, config;
                        output_prefix="plots/tmp/",
                        # outline + TM style options (surface here once)
                        show_contour::Bool=true,
                        tm_style::Union{Symbol,String}=:auto,
                        tm_len_frac::Real=0.35,
                        tm_inset_frac::Real=0.10,
                        tm_lw::Union{Real,Symbol}=:auto,
                        tm_size_frac::Real=0.28,
                        bbs = nothing)
    x = state.x_coords
    y = state.y_coords
    rot = hasproperty(state, :rotation) ? state.rotation : nothing
    box_size = state.box_size
    overlay = getfield(config, :grid_overlay)  # avoids getproperty overload surprises
    file_name = basename(config.file_path)
    N = length(x)
    monomer_path = "$(output_prefix)_$(N)_$(file_name)_placement.png"

    plot_monomers_lod(
        x, y;
        backbones = bbs,
        inputstate = state,
        rotation=rot,
        boxSize=box_size,
        monomer_radius=state.radius,     # in data units
        show_grid=overlay,
        lod=:massive,
        tick_step=:auto,
        orient_every=1,
        mode=:min,
        scale=:identity,                 # physical axes respected
        real_scale_nm=NM_PER_DATA,       # 0.37 nm/data unit
        unit=:auto,
        figsize=(1200,1200),
        dpi=200,
        save_path=monomer_path,
        draw_contacts=true,
        contact_scale=1.02,
        contact_max_lines=min(200_000, 10 * length(x)),
        show_contour=show_contour,
        tm_style=tm_style,
        tm_len_frac=tm_len_frac,
        tm_inset_frac=tm_inset_frac,
        tm_lw=tm_lw,
        tm_size_frac=tm_size_frac,
        edges = state.edges
    )

    println("Saved:\n  ", monomer_path)

    # Optional side-by-side W vs K if available on state
    if hasproperty(state, :W) && hasproperty(state, :K)
        out_wk = replace(monomer_path, "_placement.png" => "_W_vs_K.png")
        try
            plot_W_vs_K_max(state.W, state.K; save_path=out_wk)
        catch err
            @warn "plot_W_vs_K failed: $err"
        end
    end
end


############################
# Polyline Backbone Analysis
############################

# -- helpers: geometry --
hyp2(dx,dy) = sqrt(dx*dx + dy*dy)


############################
# Standalone backbone overlay saver
############################


# If you want to run directly, add your loader here:
if abspath(PROGRAM_FILE) == @__FILE__
    println("Run this via your main that provides `state, config`")
end
