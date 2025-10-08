# run_5_trials_faceCDF_and_strands.jl
# Usage:
#   julia run_5_trials_faceCDF_and_strands.jl
#
# Outputs:
#   figs/cdf/face_area_cdf_overlay.png
#   figs/strands/strand_seed_<SEED>.png  (5 images, one per trial)

using Random
using Dates
using Printf
using Statistics
using Plots

# ---------------- Project includes (adjust paths if needed) ----------------
# Keep these to the minimum needed to avoid auto-running main()
include("dependencies.jl")     # if your project has it
include("initialization.jl")   # defines Config or setup helpers
include("simulation.jl")       # run placement to get SimulationState
include("utils.jl")            # misc helpers you rely on
include("faces.jl")            # compute_enclosed_faces_xy / plot_face_area_cdf helpers
include("longest_strand.jl")   # compute_backbone
# include("plot_results.jl")   # if generate_plots is there, uncomment this

# ---------------- Config (mirror your main.jl defaults as needed) ----------
# Customize these to match your dataset + performance needs.
const FILE_PATH = "Claudins/8Cpop_grids.txt" # path to your W input
const N_MONOMERS = 50_000                     # target monomers per trial
const GRID_SIZE  = 2000
const BOX_SIZE   = 10
const NN_RESTR   = 3
const BOX_CAP    = 10
const RADIUS     = 1.0
const FT         = 5e-4                       # threshold for W floor
const GRID_OVERLAY = false

# If your Config also needs extra fields, add them here.
config = Config(
    ft=FT,
    box_size=BOX_SIZE,
    nn_restriction=NN_RESTR,
    box_capacity=BOX_CAP,
    monomer_radius=RADIUS,
    grid_size=GRID_SIZE,
    max_monomers=N_MONOMERS,
    file_path=FILE_PATH,
    grid_overlay=GRID_OVERLAY
)

# ----------------- I/O folders --------------------------------------------
mkpath("figs/cdf")
mkpath("figs/strands")
mkpath("logs")

# ----------------- Helpers: run one trial ----------------------------------
"""
    run_one_trial(seed) -> (state, areas)

Runs a single placement with RNG `seed`.
Returns the final `state` and the vector of enclosed face areas (Float64).
"""
function run_one_trial(seed::Int)
    rng = MersenneTwister(seed)
    @info "Trial with seed=$seed: starting placement…"

    # --- run your placement ---
    Random.seed!(rng)
    state = run(config)   # ✅ this is your project’s main simulation entry point

    # --- compute enclosed faces + areas ---
    areas = Float64[]
    if isdefined(@__MODULE__, :compute_enclosed_faces_xy)
        x = state.x_coords
        y = state.y_coords
        edges = hasproperty(state, :edges) ? state.edges : getproperty(state, :edge_list)
        res = compute_enclosed_faces_xy(
            x, y, edges;
            area_floor=5.0,
            drop_outer=true,
            normalize_orientation=true,
            return_abs=true
        )
        areas = res.areas
    else
        @warn "compute_enclosed_faces_xy not found; areas will be empty."
    end

    return state, areas
end


# ----------------- Helpers: CDF drawing ------------------------------------
"""
    ecdf_sorted(v) -> (xs, ys)

Return x-values = sorted(v) and y-values = empirical CDF in [0,1].
"""
function ecdf_sorted(v::AbstractVector{<:Real})
    if isempty(v)
        return Float64[], Float64[]
    end
    xs = sort!(collect(v))
    n  = length(xs)
    ys = collect(1:n) ./ n
    return xs, ys
end

"""
    plot_face_area_cdf_overlay!(plt, areas; label)

Use your project's `plot_face_area_cdf!` if it exists. Otherwise, fall back to
a simple ECDF of |area|.
"""
function plot_face_area_cdf_overlay!(plt::Plots.Plot, areas::Vector{<:Real}; label="")
    if isdefined(@__MODULE__, :plot_face_area_cdf!)
        try
            plot_face_area_cdf!(plt, areas; label=label)  # your native helper
            return
        catch err
            @warn "plot_face_area_cdf! failed, falling back to ECDF" exception=(err, catch_backtrace())
        end
    end
    xs, ys = ecdf_sorted(areas)
    plot!(plt, xs, ys; lw=2, label=label)
end

# ----------------- Backbone + Strand plotting per trial --------------------
"""
    save_strand_plot(state, seed)

Computes a backbone with `compute_backbone` and then saves a strand figure.
Tries a few signature variants for `generate_plots`.
"""
function save_strand_plot(state, seed::Int; λ::Float64=0.5, mode::Symbol=:both, use_peel::Bool=true)
    # Compute backbone
    backbone = nothing
    if isdefined(@__MODULE__, :compute_backbone)
        backbone = compute_backbone(state; λ=λ, mode=mode, use_peel=use_peel)
    else
        @warn "compute_backbone not found; saving strand without backbone overlay."
    end

    # Save using your plotter
    out = @sprintf("figs/strands/strand_seed_%d.png", seed)
    if isdefined(@__MODULE__, :generate_plots)
        ok = false
        # Try common calling styles to support your local version.
        try
            generate_plots(state; backbone=backbone, save_path=out)
            ok = true
        catch
            try
                generate_plots(state, backbone, out)
                ok = true
            catch
                try
                    generate_plots(state; save_path=out)  # no backbone path
                    ok = true
                catch err
                    @warn "generate_plots could not be called successfully." exception=(err, catch_backtrace())
                end
            end
        end
        if ok
            @info "Saved strand plot: $out"
            return
        end
    end

    # Minimal fallback: scatter monomers + (optional) backbone
    try
        x = state.x_coords; y = state.y_coords
        plt = plot(aspect_ratio=:equal, legend=:topright, grid=false, title="Strand (seed=$seed)")
        scatter!(plt, x, y; ms=1.5, alpha=0.7, label="monomers")
        if backbone !== nothing && hasproperty(state, :edges)
            # draw backbone edges as a polyline if your backbone returns path indices
            if isa(backbone, Vector{Int})
                plot!(plt, x[backbone], y[backbone]; lw=2, label="backbone")
            end
        end
        savefig(plt, out)
        @info "Saved fallback strand plot: $out"
    catch err
        @warn "Fallback strand plotting failed" exception=(err, catch_backtrace())
    end
end

# ----------------- Main: 5 trials + overlay CDF ----------------------------
function main()
    seeds = [101, 202, 303, 404, 505]  # change if you like
    all_areas = Vector{Vector{Float64}}(undef, length(seeds))
    states    = Vector{Any}(undef, length(seeds))

    # Run trials
    for (i, s) in enumerate(seeds)
        t0 = time()
        state, areas = run_one_trial(s)
        states[i] = state
        all_areas[i] = areas
        dt = time() - t0
        @info "Trial seed=$s done: faces=$(length(areas)) in $(round(dt; digits=2))s"

        # Save strand plot with backbone
        save_strand_plot(state, s)
    end

    # Build overlay plot of face-area CDFs
    plt = plot(
        xlabel="Enclosed face area",
        ylabel="CDF",
        title="Face-area CDF (5 seeds)",
        legend=:bottomright,
        grid=false
    )
    for (i, s) in enumerate(seeds)
        label = "seed=$s"
        plot_face_area_cdf_overlay!(plt, all_areas[i]; label=label)
    end

    out = "figs/cdf/face_area_cdf_overlay.png"
    savefig(plt, out)
    @info "Saved overlay CDF: $out"
end

# ----------------- Run ------------------------------------------------------
main()
