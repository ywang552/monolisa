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
function add_orientation_arcs!(
    plt, xT, yT, rotation_deg, r_plot;
    indices::AbstractVector{Int}, every::Int=50,
    arc_span_deg::Real=20.0, thickness_frac::Real=0.45,
    color=:green, alpha=0.95
)
    for j in 1:every:length(indices)
        i = indices[j]
        θ = rotation_deg[i]
        shp = ring_sector_shape(xT[i], yT[i], r_plot, thickness_frac, θ;
                                span=arc_span_deg, ns=12)
        plot!(plt, shp, c=color, opacity=alpha, linecolor=:transparent)
    end
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
    orient_len::Union{Real,Symbol}=:auto,   # (kept for fallback arrows)
    marker_size::Union{Real,Symbol}=:auto,
    figsize::Tuple{Int,Int}=(1200,1200),
    dpi::Int=200,
    # physical units
    real_scale_nm::Real = NM_PER_DATA,  # nm per data unit
    unit::Union{Symbol,String} = :auto, # :auto | :nm | :μm | :mm
    save_path=nothing
)
    N = length(x)
    if isempty(x) || isempty(y)
        @warn "No points to plot"
        return
    end

    # 1) Geometric transform (shift/fit) — ONE time
    xT, yT, s = transform_coords(x, y; mode=mode, scale=scale, target_max=target_max)

    # 2) Physical conversion to nm, then pick display unit and divide
    x_nm = xT .* real_scale_nm
    y_nm = yT .* real_scale_nm
    r_nm = monomer_radius === nothing ? nothing : Float64(monomer_radius) * s * real_scale_nm

    spanx_nm = maximum(x_nm) - minimum(x_nm)
    spany_nm = maximum(y_nm) - minimum(y_nm)
    span_nm  = max(spanx_nm, spany_nm)

    u   = unit === :auto ? pick_unit(span_nm) : Symbol(unit)
    div = unit_divisor(u)  # nm → chosen unit

    xT = x_nm ./ div
    yT = y_nm ./ div
    r_plot_unit = r_nm === nothing ? nothing : r_nm / div

    # 3) LOD
    if lod == :auto
        lod = N ≤ 100 ? :detail : (N ≤ 10_000 ? :medium : :massive)
    end

    # 4) Axis limits + ticks in chosen unit
    spanx = maximum(xT) - minimum(xT)
    spany = maximum(yT) - minimum(yT)

    if xlim === nothing || ylim === nothing
        xlim === nothing && (xlim = (minimum(xT), maximum(xT)))
        ylim === nothing && (ylim = (minimum(yT), maximum(yT)))
    end
    xticks = yticks = nothing
    if tick_step === :auto
        step = nice_tick_step(max(last(xlim)-first(xlim), last(ylim)-first(ylim)))
        x0 = first(xlim) - mod(first(xlim), step)
        y0 = first(ylim) - mod(first(ylim), step)
        xticks = collect(x0:step:last(xlim))
        yticks = collect(y0:step:last(ylim))
    elseif tick_step !== nothing
        step = Float64(tick_step)
        xticks = collect(first(xlim):step:last(xlim))
        yticks = collect(first(ylim):step:last(ylim))
    end

    # 5) Auto marker size from physical radius if available
    ms_from_radius = nothing
    if r_plot_unit !== nothing && lod != :massive
        # data radius (chosen unit) -> screen points
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


    # 6) Orientation sampling density
    target_arrows = lod == :detail ? 200 : 120
    oe_auto = max(1, Int(ceil(N / target_arrows)))
    oe = orient_every === :auto ? oe_auto : Int(orient_every)

    # 7) Plot
    plt = nothing
    if lod == :massive
        counts, xedges, yedges = density_grid(
            clamp.(xT, first(xlim), last(xlim)),
            clamp.(yT, first(ylim), last(ylim)); bins=bins
        )
        plt = heatmap(xedges, yedges, counts;
                      aspect_ratio=:equal, legend=false, dpi=dpi, size=figsize,
                      xlim=xlim, ylim=ylim, xticks=xticks, yticks=yticks,
                      color=cgrad([:blue, :red]))
    else
        idx = 1:length(xT)
        if lod == :medium && N > sample_cap_points
            stride = ceil(Int, N / sample_cap_points)
            idx = 1:stride:length(xT)
        end
        plt = scatter(xT[idx], yT[idx];
            marker=:circle, ms=ms_final,
            markercolor=:green,           # fill
            markerstrokecolor=:green,       # contour
            markerstrokewidth=0.8,
            linecolor=:transparent,
            aspect_ratio=:equal, legend=false, grid=false,
            dpi=dpi, size=figsize,
            xlim=xlim, ylim=ylim, xticks=xticks, yticks=yticks)

        # Spin ticks (rim segments) if we have rotation + radius
        if show_orientation && rotation !== nothing && !isempty(rotation)
            maxidx = maximum(collect(idx))
            if length(rotation) ≥ maxidx
                vis_ids = collect(idx)
                if r_plot_unit !== nothing
                    add_orientation_arcs!(plt, xT, yT, rotation, r_plot_unit;
                                          indices=vis_ids,
                                          every=oe,
                                          arc_span_deg=25.0,      # length along rim
                                          thickness_frac=0.45,    # thickness toward center
                                          color=:red, alpha=0.95)
                else
                    # fallback arrows if no radius available
                    oidx = vis_ids[1:oe:length(vis_ids)]
                    # small arrows, scaled to axis span
                    spanu = max(spanx, spany)
                    alen = 0.03 * spanu
                    u = cosd.(rotation[oidx]); v = sind.(rotation[oidx])
                    quiver!(plt, xT[oidx], yT[oidx], quiver=(u .* alen, v .* alen);
                            lw=0.5, alpha=0.8, linecolor=:black)
                end
            end
        end
    end

    # Optional grid overlay using same unit scale
    if show_grid && boxSize !== nothing
        # boxSize is in data units → convert to chosen unit
        scaled_box = Float64(boxSize) * s * (real_scale_nm / div)
        x_end, y_end = last(xlim), last(ylim)
        for xg in 0:scaled_box:x_end
            plot!(plt, [xg, xg], [first(ylim), y_end], lw=0.5, alpha=0.3, linecolor=:gray)
        end
        for yg in 0:scaled_box:y_end
            plot!(plt, [first(xlim), x_end], [yg, yg], lw=0.5, alpha=0.3, linecolor=:gray)
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

# -------------------------
# Public entry that uses state/config
# -------------------------
using Base.Filesystem: basename

function generate_plots(state::AbstractState, config; output_prefix="results/monolisa")
    x = state.x_coords
    y = state.y_coords
    rot = hasproperty(state, :rotation) ? state.rotation : nothing
    box_size = state.box_size
    overlay = config.grid_overlay
    file_name = basename(config.file_path)
    N = length(x)

    lod = N ≤ 100 ? :detail : (N ≤ 10_000 ? :medium : :massive)

    monomer_path = "$(output_prefix)_$(N)_$(file_name)_placement.png"
    plot_monomers_lod(
        x, y;
        rotation=rot,
        boxSize=box_size,
        monomer_radius=state.radius,   # data units
        show_grid=overlay,
        show_orientation=(lod != :massive),
        lod=lod,
        tick_step=:auto,
        orient_every = 1,
        mode=:min,
        scale=:identity,               # important: no extra fit if you want physical axes
        real_scale_nm=NM_PER_DATA,     # 0.37 nm/data unit
        unit=:auto,                    # auto-pick nm / μm / mm
        figsize=(1200,1200),
        dpi=200,
        save_path=monomer_path
    )

    println("Saved:\n  ", monomer_path)
end

# If you want to run directly, add your loader here:
if abspath(PROGRAM_FILE) == @__FILE__
    println("Run this via your main that provides `state, config`")
end
