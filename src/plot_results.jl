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
    # NEW: backbone overlay controls
    backbone_overlay::Symbol = :none,       # :none | :curvature
    overlay_corner_angle_deg::Real = 25.0,
    overlay_simplify::Union{Nothing,Real,Symbol} = :auto,  # :auto uses r_plot_unit
    overlay_resample::Union{Nothing,Real,Symbol} = :auto,  # :auto uses r_plot_unit
    overlay_lw::Real = 3.0,
    overlay_alpha::Real = 0.95
)
    N = length(x)
    if isempty(x) || isempty(y)
        @warn "No points to plot"
        return
    end

    # 1) transform once
    xT, yT, s = transform_coords(x, y; mode=mode, scale=scale, target_max=target_max)

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
        lod = N ≤ 100 ? :detail : (N ≤ 10_000 ? :medium : :massive)
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

        # compute backbone once here (we'll use it for overlay and debug PNG)
        bb_edges = nothing
        if r_plot_unit !== nothing
            bb_edges = backbone_edges_mst(xT, yT, r_plot_unit; xlim=xlim, ylim=ylim, contact_scale=contact_scale)
            # optional tiny debug image of just the backbone
            plt_backbone = plot(size=(800,800), dpi=200, aspect_ratio=:equal, legend=false)
            scatter!(plt_backbone, xT, yT; marker=:circle, ms=2.5, color=:black, label=false)
            draw_backbone!(plt_backbone, xT, yT, bb_edges; color=:red, linewidth=2.0, alpha=0.8)
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

        # TM markers (unchanged)
        if show_orientation && rotation !== nothing && !isempty(rotation)
            maxidx = maximum(collect(idx))
            if length(rotation) ≥ maxidx
                vis_ids = collect(idx)
                if r_plot_unit !== nothing
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
                    oidx = vis_ids[1:oe:length(vis_ids)]
                    spanu = max(spanx, spany)
                    alen = 0.03 * spanu
                    u = cosd.(rotation[oidx]); v = sind.(rotation[oidx])
                    quiver!(plt, xT[oidx], yT[oidx], quiver=(u .* alen, v .* alen);
                            lw=0.5, alpha=0.8, linecolor=:black, label=false)
                end
            end
        end

        # >>> NEW: curvature-colored backbone overlay <<<
        if backbone_overlay != :none && r_plot_unit !== nothing && bb_edges !== nothing
            ϵ = overlay_simplify === :auto ? r_plot_unit : (overlay_simplify === nothing ? nothing : Float64(overlay_simplify))
            Δ = overlay_resample === :auto ? r_plot_unit : (overlay_resample === nothing ? nothing : Float64(overlay_resample))
            overlay_backbone_curvature!(
                plt, xT, yT, bb_edges;
                comp_id=:all,
                simplify_eps=ϵ,
                resample_ds=Δ,
                corner_angle_deg=overlay_corner_angle_deg,
                lw=overlay_lw,
                alpha=overlay_alpha
            )
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
        backbone_overlay = :curvature,           # turn it on
        overlay_corner_angle_deg = 25.0,         # tweak if needed
        overlay_simplify = :auto,                # uses r_plot_unit
        overlay_resample = :auto                 # uses r_plot_unit

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


############################
# Polyline Backbone Analysis
############################

# -- helpers: geometry --
hyp2(dx,dy) = sqrt(dx*dx + dy*dy)

# Build adjacency and degrees on backbone nodes
function _build_adj(N::Int, edges::Vector{Tuple{Int,Int}})
    adj = [Int[] for _ in 1:N]
    for (i,j) in edges
        push!(adj[i], j); push!(adj[j], i)
    end
    return adj, [length(adj[i]) for i in 1:N]
end

# Connected components (ignore isolated nodes with deg=0)
function _components(adj)
    N = length(adj)
    seen = falses(N)
    comps = Vector{Vector{Int}}()
    for s in 1:N
        if !seen[s] && !isempty(adj[s])
            st = [s]; seen[s] = true; comp = Int[]
            while !isempty(st)
                v = pop!(st); push!(comp, v)
                for u in adj[v]
                    if !seen[u]
                        seen[u] = true; push!(st, u)
                    end
                end
            end
            push!(comps, comp)
        end
    end
    return comps
end

# Extract polylines by walking “degree-2 corridors” between terminals (deg != 2)
# Returns: Vector of polylines, each as a Vector{Int} node indices in order
function extract_polylines(N::Int, edges::Vector{Tuple{Int,Int}})
    adj, deg = _build_adj(N, edges)
    comps = _components(adj)
    visited_edge = Set{Tuple{Int,Int}}()
    polylines = Vector{Vector{Int}}()

    follow_chain = function(u::Int, v::Int)
        path = [u, v]; prev, cur = u, v
        while length(adj[cur]) == 2
            nxt = adj[cur][1] == prev ? adj[cur][2] : adj[cur][1]
            push!(path, nxt); prev, cur = cur, nxt
        end
        return path
    end

    for comp in comps
        # terminals: degree != 2 within this component
        terms = [v for v in comp if deg[v] != 2]
        # corridor walks starting at terminals
        for u in terms
            for v in adj[u]
                if !((u,v) in visited_edge || (v,u) in visited_edge)
                    path = follow_chain(u, v)
                    for k in 1:length(path)-1
                        push!(visited_edge, (path[k], path[k+1]))
                        push!(visited_edge, (path[k+1], path[k]))
                    end
                    push!(polylines, path)
                end
            end
        end
        # rare cycle fallback: if no terminals (all deg=2), walk once around
        if isempty(terms)
            v0 = first(comp); v1 = adj[v0][1]
            cyc = follow_chain(v0, v1)
            push!(polylines, cyc)
        end
    end
    return polylines
end

# Douglas–Peucker simplification on coordinates
function simplify_polyline(xs::Vector{Float64}, ys::Vector{Float64}, ε::Float64)
    n = length(xs)
    keep = falses(n); keep[1] = true; keep[n] = true
    stack = [(1,n)]
    while !isempty(stack)
        i,j = pop!(stack)
        xi, yi = xs[i], ys[i]; xj, yj = xs[j], ys[j]
        dx, dy = xj - xi, yj - yi; denom = hypot(dx,dy)
        maxd, idx = -1.0, -1
        for k in i+1:j-1
            # point-line distance
            num = abs(dy*xs[k] - dx*ys[k] + xj*yi - xi*yj)
            d = denom == 0 ? hypot(xs[k]-xi, ys[k]-yi) : num/denom
            if d > maxd; maxd = d; idx = k; end
        end
        if maxd > ε && idx > 0
            keep[idx] = true
            push!(stack, (i,idx)); push!(stack, (idx,j))
        else
            keep[i] = true; keep[j] = true
        end
    end
    return xs[keep], ys[keep]
end

# Equal-arclength resampling (chord-length parameterization)
function resample_polyline(xs::Vector{Float64}, ys::Vector{Float64}, Δs::Float64)
    n = length(xs)
    if n ≤ 2 || Δs ≤ 0; return xs, ys; end
    s = zeros(Float64, n);  # cumulative
    for k in 2:n
        s[k] = s[k-1] + hyp2(xs[k]-xs[k-1], ys[k]-ys[k-1])
    end
    L = s[end]; if L == 0; return [xs[1], xs[end]], [ys[1], ys[end]]; end
    outN = max(2, Int(floor(L/Δs))+1)
    sout = range(0, L; length=outN)
    xout = similar(sout, Float64); yout = similar(sout, Float64)
    j = 1
    for (i, si) in enumerate(sout)
        while j < n && s[j+1] < si; j += 1; end
        if j == n
            xout[i] = xs[end]; yout[i] = ys[end]
        else
            t = (si - s[j]) / max(s[j+1]-s[j], eps())
            xout[i] = (1-t)*xs[j] + t*xs[j+1]
            yout[i] = (1-t)*ys[j] + t*ys[j+1]
        end
    end
    return collect(xout), collect(yout)
end

# Turning angles (in radians) at interior vertices
function turning_angles(xs::Vector{Float64}, ys::Vector{Float64})
    n = length(xs); if n < 3; return Float64[]; end
    thetas = Float64[]
    for k in 2:n-1
        ax, ay = xs[k]-xs[k-1], ys[k]-ys[k-1]
        bx, by = xs[k+1]-xs[k], ys[k+1]-ys[k]
        na = hypot(ax,ay); nb = hypot(bx,by)
        if na == 0 || nb == 0
            push!(thetas, 0.0)
            continue
        end
        cosang = clamp((ax*bx + ay*by)/(na*nb), -1.0, 1.0)
        θ = acos(cosang)
        # signed via cross product (optional): sign = sign(ax*by - ay*bx)
        push!(thetas, θ)
    end
    return thetas
end

# Core per-strand metrics + corner detection
struct StrandMetrics
    strand_id::Int
    comp_id::Int
    n_pts::Int
    L::Float64
    D::Float64
    tortuosity::Float64
    total_abs_turn_rad::Float64
    curv_density_rad_per_unit::Float64
    mean_abs_turn_rad::Float64
    max_abs_turn_rad::Float64
    corner_count::Int
    corner_density_per_unit::Float64
    corner_s_positions::Vector{Float64}
end

function analyze_strand(xs::Vector{Float64}, ys::Vector{Float64};
                        strand_id::Int, comp_id::Int,
                        corner_angle_deg::Float64=25.0)
    n = length(xs)
    # arc length + cumulative
    s = zeros(Float64, n)
    for k in 2:n; s[k] = s[k-1] + hyp2(xs[k]-xs[k-1], ys[k]-ys[k-1]); end
    L = s[end]
    D = hyp2(xs[end]-xs[1], ys[end]-ys[1])
    τ = L > 0 ? L / max(D, eps()) : 1.0

    θs = turning_angles(xs, ys)
    absθ = abs.(θs)
    total_abs = sum(absθ)
    curv_density = L > 0 ? total_abs / L : 0.0
    mean_abs = isempty(absθ) ? 0.0 : mean(absθ)
    max_abs = isempty(absθ) ? 0.0 : maximum(absθ)

    # corner detection: angle threshold + local maxima
    thr = deg2rad(corner_angle_deg)
    corners = Int[]
    for k in 1:length(absθ)
        is_peak = (k == 1 || absθ[k] ≥ absθ[k-1]) && (k == length(absθ) || absθ[k] ≥ absθ[k+1])
        if absθ[k] ≥ thr && is_peak
            push!(corners, k+1)  # interior angle at vertex k+1
        end
    end
    corner_s = [s[i] for i in corners]
    ccount = length(corners)
    cdens = L > 0 ? ccount / L : 0.0

    return StrandMetrics(strand_id, comp_id, n, L, D, τ, total_abs, curv_density, mean_abs, max_abs, ccount, cdens, corner_s)
end

# Orchestrator: from backbone edges + (x,y) → polylines → [simplify] → [resample] → metrics
function analyze_backbone_polylines(x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                                    edges::Vector{Tuple{Int,Int}};
                                    r_plot_unit::Union{Nothing,Real}=nothing,
                                    simplify_eps::Union{Nothing,Real}=nothing,
                                    resample_ds::Union{Nothing,Real}=nothing,
                                    corner_angle_deg::Float64=25.0)
    N = length(x)
    polys = extract_polylines(N, edges)

    # map node index → component id
    adj, _ = _build_adj(N, edges)
    comps = _components(adj)
    node2comp = fill(0, N)
    for (cid, comp) in enumerate(comps)
        for v in comp; node2comp[v] = cid; end
    end

    # Prepare outputs
    strands = StrandMetrics[]
    # (Optional) component aggregates
    comp_lengths = Dict{Int,Float64}()

    for (sid, path) in enumerate(polys)
        xs = Float64[x[i] for i in path]
        ys = Float64[y[i] for i in path]
        cid = node2comp[path[1]]

        # optional simplify
        if simplify_eps !== nothing
            xs, ys = simplify_polyline(xs, ys, Float64(simplify_eps))
        end
        # optional resample
        if resample_ds !== nothing
            xs, ys = resample_polyline(xs, ys, Float64(resample_ds))
        end

        sm = analyze_strand(xs, ys; strand_id=sid, comp_id=cid, corner_angle_deg=corner_angle_deg)
        push!(strands, sm)
        comp_lengths[cid] = get(comp_lengths, cid, 0.0) + sm.L
    end

    # component summary table (simple but useful)
    comp_summary = Dict{Int,Dict{Symbol,Any}}()
    for (cid, _) in comp_lengths
        s_in = filter(sm->sm.comp_id==cid, strands)
        len_vals = [sm.L for sm in s_in]
        τ_vals   = [sm.tortuosity for sm in s_in]
        curvD    = [sm.curv_density_rad_per_unit for sm in s_in]
        comp_summary[cid] = Dict(
            :n_strands => length(s_in),
            :total_length => sum(len_vals),
            :longest_strand => maximum(len_vals),
            :median_strand => length(len_vals)>0 ? median(len_vals) : 0.0,
            :mean_tau => length(τ_vals)>0 ? mean(τ_vals) : 1.0,
            :sd_tau   => length(τ_vals)>1 ? std(τ_vals) : 0.0,
            :median_curv_density => length(curvD)>0 ? median(curvD) : 0.0,
            :mean_corners_per_unit => mean([sm.corner_density_per_unit for sm in s_in]),
        )
    end

    return strands, comp_summary
end

# (Optional) tiny overlay to inspect one component’s polylines colored by |turn angle|
# (Hook this into your plotting pipeline later if/when you want visuals.)
# using Plots
function preview_component!(plt, x, y, edges; comp_id=1, simplify_eps=nothing, resample_ds=nothing)
    N = length(x); polys = extract_polylines(N, edges)
    adj,_ = _build_adj(N, edges); comps = _components(adj)
    node2comp = fill(0, N); for (cid,comp) in enumerate(comps); for v in comp; node2comp[v]=cid; end; end
    for path in polys
        cid = node2comp[path[1]]
        cid == comp_id || continue
        xs = Float64[x[i] for i in path]; ys = Float64[y[i] for i in path]
        if simplify_eps !== nothing; xs,ys = simplify_polyline(xs,ys,simplify_eps); end
        if resample_ds !== nothing; xs,ys = resample_polyline(xs,ys,resample_ds); end
        plot!(plt, xs, ys; seriestype=:path, color=:black, alpha=0.8, lw=2)
    end
end


############################
# --- helpers for overlay ---
############################
# Build adjacency and degrees on backbone nodes
function _build_adj(N::Int, edges::Vector{Tuple{Int,Int}})
    adj = [Int[] for _ in 1:N]
    for (i,j) in edges
        push!(adj[i], j); push!(adj[j], i)
    end
    return adj, [length(adj[i]) for i in 1:N]
end

# Connected components (ignore isolated nodes with deg=0)
function _components(adj::Vector{Vector{Int}})
    N = length(adj)
    seen = falses(N)
    comps = Vector{Vector{Int}}()
    for s in 1:N
        if !seen[s] && !isempty(adj[s])
            st = [s]; seen[s] = true; comp = Int[]
            while !isempty(st)
                v = pop!(st); push!(comp, v)
                for u in adj[v]
                    if !seen[u]
                        seen[u] = true; push!(st, u)
                    end
                end
            end
            push!(comps, comp)
        end
    end
    return comps
end

# Extract polylines by walking degree-2 corridors between terminals (deg != 2)
function extract_polylines(N::Int, edges::Vector{Tuple{Int,Int}})
    adj, deg = _build_adj(N, edges)
    comps = _components(adj)
    visited_edge = Set{Tuple{Int,Int}}()
    polylines = Vector{Vector{Int}}()

    follow_chain = function(u::Int, v::Int)
        path = [u, v]; prev, cur = u, v
        while length(adj[cur]) == 2
            nxt = adj[cur][1] == prev ? adj[cur][2] : adj[cur][1]
            push!(path, nxt); prev, cur = cur, nxt
        end
        return path
    end

    for comp in comps
        terms = [v for v in comp if deg[v] != 2]
        for u in terms
            for v in adj[u]
                if !((u,v) in visited_edge || (v,u) in visited_edge)
                    path = follow_chain(u, v)
                    for k in 1:length(path)-1
                        push!(visited_edge, (path[k], path[k+1]))
                        push!(visited_edge, (path[k+1], path[k]))
                    end
                    push!(polylines, path)
                end
            end
        end
        if isempty(terms) && !isempty(comp)  # rare cycles
            v0 = first(comp); v1 = adj[v0][1]
            push!(polylines, follow_chain(v0, v1))
        end
    end
    return polylines
end

# Turning angles (radians) at interior vertices of a polyline
function turning_angles(xs::Vector{Float64}, ys::Vector{Float64})
    n = length(xs); if n < 3; return Float64[]; end
    thetas = Float64[]
    for k in 2:n-1
        ax, ay = xs[k]-xs[k-1], ys[k]-ys[k-1]
        bx, by = xs[k+1]-xs[k], ys[k+1]-ys[k]
        na = hypot(ax,ay); nb = hypot(bx,by)
        if na == 0 || nb == 0
            push!(thetas, 0.0); continue
        end
        cosang = clamp((ax*bx + ay*by)/(na*nb), -1.0, 1.0)
        push!(thetas, acos(cosang))   # unsigned; use sign via cross if needed
    end
    return thetas
end

# (re)use your simplify/resample if you already added them; otherwise light no-ops:
if !@isdefined simplify_polyline
    simplify_polyline(xs::Vector{Float64}, ys::Vector{Float64}, ε::Float64) = (xs, ys)
end
if !@isdefined resample_polyline
    resample_polyline(xs::Vector{Float64}, ys::Vector{Float64}, Δs::Float64) = (xs, ys)
end

# map interior angles to per-segment values for coloring
function _segment_turn_values(xs::Vector{Float64}, ys::Vector{Float64})
    n = length(xs)
    if n < 3
        return Float64[], Float64[]
    end
    θ = turning_angles(xs, ys)          # length n-2, at vertices 2..n-1
    absθ = abs.(θ)
    segv = Vector{Float64}(undef, n-1)
    segv[1] = absθ[1]
    for k in 2:n-2
        segv[k] = absθ[k-1]
    end
    segv[n-1] = absθ[end]
    return absθ, segv
end

# corner picker: threshold + local max in absθ (angles live at vertices 2..n-1)
function _corner_indices_from_absθ(absθ::Vector{Float64}, thr_rad::Float64)
    corners = Int[]
    m = length(absθ)
    for i in 1:m
        is_peak = (i == 1 || absθ[i] ≥ absθ[i-1]) && (i == m || absθ[i] ≥ absθ[i+1])
        if absθ[i] ≥ thr_rad && is_peak
            push!(corners, i + 1)  # map to vertex index
        end
    end
    return corners
end

# Curvature-colored overlay with corner markers
function overlay_backbone_curvature!(
    plt,
    x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
    edges::Vector{Tuple{Int,Int}};
    comp_id::Union{Int,Symbol} = :all,
    simplify_eps::Union{Nothing,Real} = nothing,
    resample_ds::Union{Nothing,Real} = nothing,
    corner_angle_deg::Real = 25.0,
    # visuals
    cmap_name = :viridis,
    lw::Real = 3.0,
    alpha::Real = 0.95,
    corner_marker = :xcross,
    corner_ms_min::Real = 4.0,
    corner_ms_max::Real = 12.0,
    corner_color = :red,
    show_colorbar::Bool = true,
    show_component_labels::Bool = true,
    component_label_fontsize::Int = 9,
)
    N = length(x)
    polys = extract_polylines(N, edges)

    # build components + map node -> component id
    adj, _ = _build_adj(N, edges)
    comps = _components(adj)
    node2comp = fill(0, N)
    for (cid, comp) in enumerate(comps), v in comp
        node2comp[v] = cid
    end

    # pass 1: prepare polylines and collect segment curvature (in degrees)
    all_seg_vals_deg = Float64[]
    # cache entries: (comp_id, xs, ys, absθ_deg, segv_deg)
    cache = Vector{Tuple{Int,Vector{Float64},Vector{Float64},Vector{Float64},Vector{Float64}}}()
    for path in polys
        cid = node2comp[path[1]]
        if comp_id !== :all && cid != comp_id
            continue
        end

        xs = Float64[x[i] for i in path]
        ys = Float64[y[i] for i in path]
        if simplify_eps !== nothing
            xs, ys = simplify_polyline(xs, ys, Float64(simplify_eps))
        end
        if resample_ds !== nothing
            xs, ys = resample_polyline(xs, ys, Float64(resample_ds))
        end

        absθ_rad, segv_rad = _segment_turn_values(xs, ys)   # both in radians
        absθ_deg = rad2deg.(absθ_rad)
        segv_deg = rad2deg.(segv_rad)

        push!(cache, (cid, xs, ys, absθ_deg, segv_deg))
        append!(all_seg_vals_deg, segv_deg)
    end

    # color scale based on curvature in degrees
    vmin_deg, vmax_deg = isempty(all_seg_vals_deg) ? (0.0, 0.0) : extrema(all_seg_vals_deg)
    Δ = max(vmax_deg - vmin_deg, eps())
    cmap = cgrad(cmap_name)

    # pass 2: draw colored segments + corner markers
    thr_deg = corner_angle_deg
    for (_, xs, ys, absθ_deg, segv_deg) in cache
        Lseg = length(segv_deg) # number of edges in the polyline
        for k in 1:Lseg
            v = (segv_deg[k] - vmin_deg) / Δ
            col = cmap[clamp(v, 0.0, 1.0)]
            plot!(plt, xs[k:k+1], ys[k:k+1];
                  seriestype = :path, color = col, lw = lw, alpha = alpha, label = false)
        end

        # corner markers (vertex-based). absθ_deg has length = n_points - 2
        if !isempty(absθ_deg)
            # IMPORTANT: _corner_indices_from_absθ expects radians -> convert
            corner_idx = _detect_corners_multiscale(xs, ys,
                                        deg2rad.(segv_deg),  # we have deg; convert back to rad
                                        deg2rad.(absθ_deg);
                                        thr_deg = corner_angle_deg,
                                        win = 5,
                                        min_sep_pts = 3,
                                        min_arc_sep = 0.75 * mean(diff(xs)))
            for vidx in corner_idx
                # interior vertex at 2..length(xs)-1 maps to absθ index (vidx-1)
                aidx = vidx - 1
                aidx < 1 && continue
                aidx > length(absθ_deg) && continue
                frac = clamp(absθ_deg[aidx] / 180.0, 0.0, 1.0)
                msz = corner_ms_min + (corner_ms_max - corner_ms_min) * frac
                scatter!(plt, [xs[vidx]], [ys[vidx]];
                         marker = corner_marker, ms = msz, color = corner_color,
                         alpha = 1.0, label = false)
            end
        end
    end

    # optional colorbar: use a 1x2 transparent heatmap to carry the scale
    if show_colorbar
        dummy = [vmin_deg vmin_deg; vmax_deg vmax_deg]
        heatmap!(plt, [0, 1], [0, 1], dummy;
                 alpha = 0.0, color = cmap, clims = (vmin_deg, vmax_deg),
                 colorbar = true, label = false)
    end

    # optional component labels at centroids
    if show_component_labels
        sel = comp_id === :all ? (1:length(comps)) : (comp_id:comp_id)
        for cid in sel
            comp = comps[cid]
            isempty(comp) && continue
            cx = mean(@view x[comp]); cy = mean(@view y[comp])
            annotate!(plt, cx, cy, text("C$cid", component_label_fontsize, :black))
        end
    end

    return plt
end

############################
# Standalone backbone overlay saver
############################
using Plots

function save_backbone_overlay(
    state;
    # file
    save_path::AbstractString,
    figsize::Tuple{Int,Int}=(1200,1200),
    dpi::Int=200,
    # geometry / units (match plot_monomers_lod defaults)
    mode::Symbol=:min, scale::Union{Symbol,Real}=:identity, target_max=nothing,
    real_scale_nm::Real = NM_PER_DATA,
    unit::Union{Symbol,String} = :auto,
    # backbone + overlay knobs
    contact_scale::Real = 1.02,
    overlay_corner_angle_deg::Real = 25.0,
    overlay_simplify::Union{Nothing,Real,Symbol} = :auto,  # :auto uses r_plot_unit
    overlay_resample::Union{Nothing,Real,Symbol} = :auto,  # :auto uses r_plot_unit
    overlay_lw::Real = 4.0,
    overlay_alpha::Real = 1.0,
)
    x = state.x_coords
    y = state.y_coords
    monomer_radius = state.radius

    isempty(x) && isempty(y) && error("No coordinates in state")

    # --- transform + units (same as plot_monomers_lod) ---
    xT0, yT0, s = transform_coords(x, y; mode=mode, scale=scale, target_max=target_max)
    x_nm = xT0 .* real_scale_nm
    y_nm = yT0 .* real_scale_nm
    r_nm = monomer_radius === nothing ? nothing : Float64(monomer_radius) * s * real_scale_nm

    spanx_nm = maximum(x_nm) - minimum(x_nm)
    spany_nm = maximum(y_nm) - minimum(y_nm)
    span_nm  = max(spanx_nm, spany_nm)
    u   = unit === :auto ? pick_unit(span_nm) : Symbol(unit)
    div = unit_divisor(u)

    xT = x_nm ./ div
    yT = y_nm ./ div
    r_plot_unit = r_nm === nothing ? nothing : r_nm / div

    # --- limits with a little padding (like the main plot) ---
    spanx = maximum(xT) - minimum(xT)
    spany = maximum(yT) - minimum(yT)
    pad = r_plot_unit !== nothing ? 1.15 * r_plot_unit : 0.02 * max(spanx, spany)
    xlim = (minimum(xT) - pad, maximum(xT) + pad)
    ylim = (minimum(yT) - pad, maximum(yT) + pad)

    # --- compute backbone on transformed coords ---
    r_plot_unit === nothing && error("state.radius must be set for backbone overlay")
    bb_edges = backbone_edges_mst(xT, yT, r_plot_unit; xlim=xlim, ylim=ylim, contact_scale=contact_scale)

    # --- make a clean figure and draw only the curvature overlay ---
    plt = plot(; size=figsize, dpi=dpi, aspect_ratio=:equal, legend=false,
                xlim=xlim, ylim=ylim, background_color=:white)

    ϵ = overlay_simplify === :auto ? r_plot_unit : (overlay_simplify === nothing ? nothing : Float64(overlay_simplify))
    Δ = overlay_resample === :auto ? r_plot_unit : (overlay_resample === nothing ? nothing : Float64(overlay_resample))

    overlay_backbone_curvature!(
        plt, xT, yT, bb_edges;
        comp_id=:all,
        simplify_eps=.1,
        resample_ds=.3,
        corner_angle_deg=25,
        lw=4,
        alpha=overlay_alpha,
        show_colorbar = true,
        show_component_labels = false   # suppress C1, C2, ...
           
    )

    xlabel!(plt, "Distance ($(unit_label(u)))")
    ylabel!(plt, "Distance ($(unit_label(u)))")

    # save
    mkpath(dirname(save_path))
    savefig(plt, save_path)
    println("Saved backbone overlay → $save_path")
    return plt
end


"""
Robust corner detection on a resampled polyline.

Inputs:
- xs, ys: resampled polyline (length M)
- seg_turn_rad: per-segment turning (length M-1), from `_segment_turn_values`
- absθ_rad: per-vertex absolute interior angles (length M-2), from `_segment_turn_values`
- thr_deg: primary threshold on *smoothed* turning angle (degrees)
- win: smoothing window (odd, >=3) for moving average on segment turns
- min_sep_pts: non-maximum suppression window (corners must be >= this many points apart)
- min_arc_sep: minimum arc-length separation between accepted corners (in plot units)

Returns:
- corner vertex indices in the polyline indexing (2..M-1)
"""
function _detect_corners_multiscale(xs::Vector{Float64}, ys::Vector{Float64},
                                    seg_turn_rad::Vector{Float64},
                                    absθ_rad::Vector{Float64};
                                    thr_deg::Real = 25.0,
                                    win::Int = 5,
                                    min_sep_pts::Int = 3,
                                    min_arc_sep::Real = 0.75)

    M = length(xs)
    if M < 3
        return Int[]
    end

    # 1) Smooth per-segment turning (Savitzky would be nicer; MA is fine & fast)
    L = length(seg_turn_rad)
    w = max(3, isodd(win) ? win : win + 1)
    half = (w - 1) ÷ 2
    pad = [seg_turn_rad[1:half]; seg_turn_rad; seg_turn_rad[end-half+1:end]]
    sm = similar(seg_turn_rad)
    for i in 1:L
        sm[i] = mean(@view pad[i:i+w-1])
    end

    # 2) Threshold on smoothed segment turning (in degrees)
    sm_deg = rad2deg.(sm)
    mask = sm_deg .>= thr_deg

    # 3) Non-maximum suppression on the mask
    #    Keep local peaks; discard neighbors within +/- min_sep_pts if smaller
    peaks = falses(L)
    i = 1
    while i <= L
        if !mask[i]
            i += 1
            continue
        end
        # look at a local group of consecutive "true" positions
        j = i
        while j <= L && mask[j]
            j += 1
        end
        # in [i, j-1], pick the index with max smoothed value
        seg_slice = sm_deg[i:j-1]
        rel_max = argmax(seg_slice)
        k = i + rel_max - 1
        peaks[k] = true
        i = j
    end

    # spread-based NMS: remove peaks closer than min_sep_pts (by index)
    peak_idx = findall(peaks)
    filtered = Int[]
    last_keep = -10^9
    for p in peak_idx
        if p - last_keep >= min_sep_pts
            push!(filtered, p)
            last_keep = p
        end
    end

    # 4) Map per-segment peak (k in 1..M-1) to a *vertex* index in 2..M-1
    #    A corner lives at vertex v ~ k+1 (since segment k goes v..v+1)
    cand_vertices = clamp.(filtered .+ 1, 2, M-1)

    # 5) Enforce minimum arc-length separation in actual distance
    # precompute cumulative arclength
    s = zeros(Float64, M)
    for t in 2:M
        s[t] = s[t-1] + hypot(xs[t]-xs[t-1], ys[t]-ys[t-1])
    end

    kept = Int[]
    last_s = -1e18
    for v in cand_vertices
        if (s[v] - last_s) >= min_arc_sep
            push!(kept, v)
            last_s = s[v]
        end
    end

    return kept
end



# If you want to run directly, add your loader here:
if abspath(PROGRAM_FILE) == @__FILE__
    println("Run this via your main that provides `state, config`")
end
