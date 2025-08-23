############################
# Monomer Plotting (units + spin tick)
############################

using Plots
using Plots: Shape
using StatsBase

# -------------------------
# Config: physical mapping
# -------------------------
const NM_PER_DATA = 0.37  # 1 data unit = 0.37 nm

# -------------------------
# Small utilities
# -------------------------
ceil_to(x, step) = step * ceil(Int, x / step)
build_ticks(limit, step) = collect(0:step:limit)
const DEG2RAD = π/180



"""
Draw contact lines between disks of radius r at positions (x,y).
Uses a spatial grid so we only check local neighbors.
- Only draws if r is provided (not `nothing`)
- contact_scale lets you be a bit lenient (e.g., 1.02).
"""
function _draw_contacts_grid!(plt, x::AbstractVector, y::AbstractVector, r::Real;
                              xlim::Tuple, ylim::Tuple,
                              contact_scale::Real=1.02,
                              max_lines::Int=200_000)

    N = length(x)
    contact_r = 2*r*contact_scale          # <- ensure the * is here (not `2r`)
    cell = contact_r

    # grid bucketing
    gx(i) = Int(floor((x[i] - first(xlim)) / cell))
    gy(i) = Int(floor((y[i] - first(ylim)) / cell))
    buckets = Dict{Tuple{Int,Int}, Vector{Int}}()
    for i in 1:N
        key = (gx(i), gy(i))
        push!(get!(buckets, key, Int[]), i)
    end

    offs = ((-1,-1),(-1,0),(-1,1),(0,-1),(0,0),(0,1),(1,-1),(1,0),(1,1))
    drawn = 0
    csq = contact_r^2

    for i in 1:N
        xi = x[i]; yi = y[i]
        isfinite(xi) && isfinite(yi) || continue
        gi, gj = gx(i), gy(i)

        for (di, dj) in offs
            nbrs = get(buckets, (gi+di, gj+dj), nothing)
            isnothing(nbrs) && continue
            for j in nbrs
                j <= i && continue
                xj = x[j]; yj = y[j]
                isfinite(xj) && isfinite(yj) || continue

                dx = xj - xi; dy = yj - yi
                if dx*dx + dy*dy <= csq
                    # use 2-point VECTORS, not scalars/tuples
                    plot!(plt, [xi, xj], [yi, yj];
                          seriestype=:path, color=:grey, alpha=1, linewidth=4)
                    drawn += 1
                    drawn >= max_lines && return
                end
            end
        end
    end
end

# -------------------------
# Contact backbone (return-only)
# -------------------------

# Return contact edges (i,j) with their Euclidean distances (no plotting).
function find_contact_edges(x::AbstractVector, y::AbstractVector, r::Real;
                            xlim::Tuple, ylim::Tuple, contact_scale::Real=1.02)
    N = length(x)
    contact_r = 2r * contact_scale
    cell = contact_r
    csq = contact_r^2

    gx(i) = Int(floor((x[i] - first(xlim)) / cell))
    gy(i) = Int(floor((y[i] - first(ylim)) / cell))

    buckets = Dict{Tuple{Int,Int}, Vector{Int}}()
    for i in 1:N
        xi, yi = x[i], y[i]
        (isfinite(xi) && isfinite(yi)) || continue
        push!(get!(buckets, (gx(i), gy(i)), Int[]), i)
    end

    offs  = ((-1,-1),(-1,0),(-1,1),(0,-1),(0,0),(0,1),(1,-1),(1,0),(1,1))
    edges = Vector{Tuple{Int,Int}}()
    dists = Float64[]

    @inbounds for i in 1:N
        xi = x[i]; yi = y[i]
        (isfinite(xi) && isfinite(yi)) || continue
        gi, gj = gx(i), gy(i)

        for (di, dj) in offs
            nbrs = get(buckets, (gi+di, gj+dj), nothing)
            isnothing(nbrs) && continue
            for j in nbrs
                j <= i && continue
                xj = x[j]; yj = y[j]
                (isfinite(xj) && isfinite(yj)) || continue
                dx = xj - xi; dy = yj - yi
                dsq = dx*dx + dy*dy
                if dsq <= csq
                    push!(edges, (i, j))
                    push!(dists, sqrt(dsq))
                end
            end
        end
    end

    return edges, dists
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

# Optional: draw a set of edges on a given plot (overlay helper)
function draw_backbone!(plt, x::AbstractVector, y::AbstractVector, edges::Vector{Tuple{Int,Int}};
                        color=:black, alpha=0.9, linewidth=3.0)
    @inbounds for (i, j) in edges
        plot!(plt, [x[i], x[j]], [y[i], y[j]];label = false,
              seriestype=:path, color=color, alpha=alpha, linewidth=linewidth, linecap=:round)
    end
    return plt
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
Make a filled ring-segment polygon centered at (cx,cy).
- r: outer radius (plot units, i.e., final chosen unit)
- t: thickness as fraction of r (0..1)
- θ: center angle (deg; 0° = +x axis, CCW)
"""
function ring_sector_shape(cx, cy, r, t, θ; span=20.0, ns=12)
    rin  = max(r*(1 - t), 0.0)
    θ0   = (θ - span/2) * DEG2RAD
    θ1   = (θ + span/2) * DEG2RAD
    outer = [(cx + r*cos(ang),  cy + r*sin(ang)) for ang in range(θ0, θ1, length=ns)]
    inner = [(cx + rin*cos(ang), cy + rin*sin(ang)) for ang in range(θ1, θ0, length=ns)]
    pts = vcat(outer, inner)
    xs = first.(pts); ys = last.(pts)
    return Shape(xs, ys)
end

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
    draw_contacts::Union{Bool,Symbol}=:auto,      # unchanged API
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
    tm_style::Union{Symbol,String}=:auto,   # :auto | :arcs | :ticks | :triangles
    tm_len_frac::Real=0.35,                 # ticks-only
    tm_inset_frac::Real=0.10,               # ticks-only
    tm_lw::Union{Real,Symbol}=:auto,        # ticks-only
    tm_size_frac::Real=0.28                 # triangles-only
)
    N = length(x)
    if isempty(x) || isempty(y)
        @warn "No points to plot"
        return
    end

    # 1) transform once
    xT, yT, s = transform_coords(x, y; mode=mode, scale=scale, target_max=target_max)

    # 2) convert to physical units and then chosen display unit
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
        lod = N ≤ 100 ? :detail : (N ≤ 10_000 ? :medium : :massive)
    end

    # 4) limits + ticks (pad to avoid clipping)
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

    if lod == :massive
        counts, xedges, yedges = density_grid(
            clamp.(xT, first(xlim), last(xlim)),
            clamp.(yT, first(ylim), last(ylim)); bins=bins
        )
        plt = heatmap(xedges, yedges, counts;
            aspect_ratio=:equal, legend=false, dpi=dpi, size=figsize,
            xlim=xlim, ylim=ylim, xticks=xtick_vals, yticks=ytick_vals,
            color=cgrad([:blue, :red]))
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

        # optional contact lines (uses your grid accel)
        do_contacts = (draw_contacts === :auto ? (lod != :massive) : draw_contacts) && (r_plot_unit !== nothing)
        if do_contacts
            _draw_contacts_grid!(plt, xT[idx], yT[idx], r_plot_unit;
                                 xlim=xlim, ylim=ylim,
                                 contact_scale=contact_scale,
                                 max_lines=contact_max_lines)
        end

        # compass legend in least-occupied corner
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

        # --- TM marker rendering (style switch) ---
        if show_orientation && rotation !== nothing && !isempty(rotation)
            maxidx = maximum(collect(idx))
            if length(rotation) ≥ maxidx
                vis_ids = collect(idx)
                if r_plot_unit !== nothing
                    bb_edges = backbone_edges_mst(xT, yT, r_plot_unit; xlim=xlim, ylim=ylim, contact_scale=contact_scale)
                    plt_backbone = plot(size=(1200,1200), dpi=200, aspect_ratio=:equal, legend=false)
                    scatter!(plt_backbone, xT, yT; marker=:circle, ms=3, color=:black, label=false)
                    draw_backbone!(plt_backbone, xT, yT, bb_edges; color=:red, linewidth=2.0, alpha=0.8)
                    if save_path !== nothing
                        savefig(plt_backbone, replace(save_path, ".png" => "_backbone.png"))
                    end
                    style = tm_style == :auto ? (show_contour ? :arcs : :ticks) : Symbol(tm_style)

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
                    # fallback arrows if no radius
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



using Plots

"""
Plot W (probabilities) and K (counts) with different color mappings.
- W: clamped/normalized to [0,1], perceptual palette (:viridis).
- K: log10(count+1) to compress heavy tails, contrasty palette (:magma).
- Angle ticks at 0°, 90°, 180°, 270° assuming 72×72 (5° bins).
"""
function plot_W_vs_K(W::AbstractMatrix, K::AbstractMatrix;
                     normalize_W::Bool=true,
                     add_angle_ticks::Bool=true,
                     save_path::Union{Nothing,String}=nothing)

    @assert size(W) == size(K) "W and K must have same shape (e.g., 72×72)."
    ZED = cgrad(:roma, 10, categorical = true, scale = :exp)

    # --- W: probability map ---
    Wp = normalize_W ? clamp.(W, 0, 1) : W
    clW = normalize_W ? (0.0, 1.0) : (minimum(Wp), maximum(Wp))
    pltW = heatmap(Wp;
        aspect_ratio = :equal,
        color        = ZED,  # good for probabilities
        clims        = clW,
        colorbar_title = normalize_W ? "prob (0–1)" : "prob",
        title       = "Input connection probability (W)")

    # --- K: count map with log compression ---
    Klog = log10.(K .+ 1)  # log10 for readable dynamic range
    kmin, kmax = extrema(Klog)
    kmax == kmin && (kmax += 1e-6)  # avoid degenerate color scale
    pltK = heatmap(Klog;
        aspect_ratio   = :equal,
        color          = ZED, # distinct from W
        clims          = (kmin, kmax),
        colorbar_title = "log10(count + 1)",
        title          = "Simulated connections (K)")

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
using Base.Filesystem: basename

function generate_plots(state::AbstractState, config;
                        output_prefix="plots/tmp/",
                        # outline + TM style options (surface here once)
                        show_contour::Bool=true,
                        tm_style::Union{Symbol,String}=:auto,
                        tm_len_frac::Real=0.35,
                        tm_inset_frac::Real=0.10,
                        tm_lw::Union{Real,Symbol}=:auto,
                        tm_size_frac::Real=0.28)
    x = state.x_coords
    y = state.y_coords
    rot = hasproperty(state, :rotation) ? state.rotation : nothing
    box_size = state.box_size
    overlay = getfield(config, :grid_overlay)  # avoids getproperty overload surprises
    file_name = basename(config.file_path)
    N = length(x)

    lod = N ≤ 100 ? :detail : (N ≤ 10_000 ? :medium : :massive)

    monomer_path = "$(output_prefix)_$(N)_$(file_name)_placement.png"

    plot_monomers_lod(
        x, y;
        rotation=rot,
        boxSize=box_size,
        monomer_radius=state.radius,     # in data units
        show_grid=overlay,
        show_orientation=(lod != :massive),
        lod=lod,
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
        # NEW
        show_contour=show_contour,
        tm_style=tm_style,
        tm_len_frac=tm_len_frac,
        tm_inset_frac=tm_inset_frac,
        tm_lw=tm_lw,
        tm_size_frac=tm_size_frac,
    )

    println("Saved:\n  ", monomer_path)

    # Optional side-by-side W vs K if available on state
    if hasproperty(state, :W) && hasproperty(state, :K)
        out_wk = replace(monomer_path, "_placement.png" => "_W_vs_K.png")
        try
            plot_W_vs_K(state.W, state.K; normalize_W=true, save_path=out_wk)
        catch err
            @warn "plot_W_vs_K failed: $err"
        end
    end
end



# If you want to run directly, add your loader here:
if abspath(PROGRAM_FILE) == @__FILE__
    println("Run this via your main that provides `state, config`")
end
