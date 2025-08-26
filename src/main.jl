# ==============================================
# File: main.jl
# Version: v1.0.0
# Description: Entry point for MonoLisa simulation.
# Author: [Youchuan Wang]
# Date: [11/17/2024]
# ==============================================

include("dependencies.jl")
include("initialization.jl")
include("simulation.jl")
include("utils.jl")
include("plot_results.jl")


"""
Sets up the necessary folders for the project.
Creates folders if they do not already exist.
"""
function setup_project_folders()
    folders = ["saved_states", "saved_states\\minimal_states", "saved_states\\full_states", "Claudins", "plots", "logs"]
    for folder in folders
        folder = pwd()*"\\"*folder
        if !isdir(folder)
            mkdir(folder)
            println("Created folder: $folder")
        else
            println("Folder already exists: $folder")
        end
    end
end

# setup_project_folders()

config = Config(
    ft=0.0005,
    box_size=10,
    nn_restriction=3,
    box_capacity=16,
    monomer_radius=1,
    grid_size=80,
    # max_monomers=1000000,
    max_monomers=100,
    # file_path=ARGS[1],
    file_path="Claudins/C2WT_masked90_zeroint.txt",
    grid_overlay = false,
    prog = 1
)

function main()
    println("starting the placement...")
    # Pass configuration to simulation
    state = run(config)

    # Save results or generate plots
    return state
end

state = main()
print(length(state.x_coords))
out_dir = joinpath(pwd(), "plots", "tmp")
mkpath(out_dir)

name_noext, _ = splitext(basename(config.file_path))
prefix = joinpath(out_dir, name_noext*"_123")

generate_plots(state, config; output_prefix = prefix*"arcs", show_contour=true, tm_style=:arcs)
overlay_path = prefix * "overlay_curvature.png"
save_backbone_overlay(state; save_path = overlay_path,
    overlay_corner_angle_deg = 25.0,   # tweak if you want stricter/looser corners
    overlay_simplify = :auto,          # uses state.radius in plot units
    overlay_resample = :auto,          # uses state.radius in plot units
    overlay_lw = 5.0,                  # a bit thicker for visibility
    overlay_alpha = 1.0
)


# generate_plots(state, config; output_prefix = prefix*"arcs", show_contour=false, tm_style=:arcs)
# generate_plots(state, config; output_prefix = prefix*"ticks", show_contour=false, tm_style=:ticks)
# generate_plots(state, config; output_prefix = prefix*"triangles", show_contour=false, tm_style=:triangles, tm_size_frac=0.32)

# state.rotation

# x = state.x_coords
# y = state.y_coords;








# Pull coordinates + radius from state
x = state.x_coords
y = state.y_coords
r_plot_unit = state.radius   # already in plot units

# Build backbone edges first (same as you did for plotting)
bb_edges = backbone_edges_mst(
    x, y, r_plot_unit;
    xlim=(minimum(x), maximum(x)),
    ylim=(minimum(y), maximum(y)),
    contact_scale=1.02
)

# Set simplify / resample defaults relative to radius
ϵ = 1.0 * r_plot_unit     # Douglas–Peucker tolerance
Δ = 1.0 * r_plot_unit     # resample spacing

# Run analysis
strands, comps = analyze_backbone_polylines(
    x, y, bb_edges;
    r_plot_unit   = r_plot_unit,
    simplify_eps  = ϵ,
    resample_ds   = Δ,
    corner_angle_deg = 25.0
)

println("Per-strand metrics:")
for sm in strands
    println("Strand $(sm.strand_id) (comp $(sm.comp_id)): length=$(sm.L), τ=$(sm.tortuosity), corners=$(sm.corner_count)")
end

println("\nPer-component summaries:")
for (cid, summary) in comps
    println("Component $cid => ", summary)
end


# === Build tables, save, and plot summaries ===
using DataFrames, Statistics, StatsBase, CSV, JSON3, Dates, Plots

# 1) Convert your outputs into tidy DataFrames
function build_backbone_tables(strands, comps)
    # strands is a Vector of structs: (strand_id, comp_id, L, tortuosity, corner_count)
    strands_df = DataFrame(
        strand_id     = Int[],
        comp_id       = Int[],
        length        = Float64[],
        tau           = Float64[],
        corners       = Int[],
        curv_density  = Float64[],   # corners per unit length
    )
    for s in strands
        len = s.L
        c   = s.corner_count
        push!(strands_df, (
            s.strand_id,
            s.comp_id,
            len,
            s.tortuosity,
            c,
            c / max(len, eps())
        ))
    end

    # comps is Dict{Int, Dict{Symbol,Any}} (as printed above)
    components_df = DataFrame(
        component_id           = Int[],
        n_strands              = Int[],
        longest_strand         = Float64[],
        median_strand          = Float64[],
        total_length           = Float64[],
        mean_tau               = Float64[],
        sd_tau                 = Float64[],
        mean_corners_per_unit  = Float64[],
        median_curv_density    = Float64[],
    )
    for (cid, d) in comps
        push!(components_df, (
            cid,
            d[:n_strands],
            d[:longest_strand],
            d[:median_strand],
            d[:total_length],
            d[:mean_tau],
            d[:sd_tau],
            d[:mean_corners_per_unit],
            d[:median_curv_density],
        ))
    end
    return strands_df, components_df
end

# 2) Save as CSV and JSON (prefix controls filenames)
function save_backbone_metrics(strands_df::DataFrame, components_df::DataFrame;
                               outdir::AbstractString="logs/metrics",
                               prefix::AbstractString="backbone")
    isdir(outdir) || mkpath(outdir)
    ts = Dates.format(now(), "yyyymmdd_HHMM")
    base = joinpath(outdir, "$(prefix)_$(ts)")

    CSV.write(base * "_strands.csv", strands_df)
    CSV.write(base * "_components.csv", components_df)

    blob = Dict(
        "strands"    => Tables.columntable(strands_df),
        "components" => Tables.columntable(components_df),
    )
    open(base * ".json", "w") do io
        JSON3.write(io, blob; indent=2)
    end
    return (strands_csv = base * "_strands.csv",
            components_csv = base * "_components.csv",
            json = base * ".json",
            base = base)
end

# 3) Plots: histograms for strand metrics, bar charts for component metrics
function plot_backbone_metrics(strands_df::DataFrame, components_df::DataFrame;
                               outbase::AbstractString=nothing)
    # Strand-level histograms
    plt1 = histogram(strands_df.length; bins=20,
                     xlabel="Strand length", ylabel="Count",
                     title="Distribution of strand lengths", legend=false)
    plt2 = histogram(strands_df.tau; bins=20,
                     xlabel="τ (tortuosity)", ylabel="Count",
                     title="Distribution of τ", legend=false)
    plt3 = histogram(strands_df.curv_density; bins=20,
                     xlabel="Curvature density (corners/length)", ylabel="Count",
                     title="Distribution of curvature density", legend=false)

    # Strand-level scatter (optional but useful)
    plt4 = scatter(strands_df.length, strands_df.tau;
                   xlabel="Strand length", ylabel="τ",
                   title="Length vs. tortuosity", legend=false)
    plt5 = scatter(strands_df.length, strands_df.curv_density;
                   xlabel="Strand length", ylabel="Curvature density",
                   title="Length vs. curvature density", legend=false)

    # Component-level bars
    # Ensure sorted by component id for consistent x ordering
    comps_sorted = sort(components_df, :component_id)
    plt6 = bar(comps_sorted.component_id, comps_sorted.n_strands;
               xlabel="Component", ylabel="Number of strands",
               title="Strands per component", legend=false)
    plt7 = bar(comps_sorted.component_id, comps_sorted.total_length;
               xlabel="Component", ylabel="Total length",
               title="Total length per component", legend=false)
    plt8 = bar(comps_sorted.component_id, comps_sorted.mean_tau;
               xlabel="Component", ylabel="Mean τ",
               title="Mean tortuosity per component", legend=false)
    plt9 = bar(comps_sorted.component_id, comps_sorted.mean_corners_per_unit;
               xlabel="Component", ylabel="Mean corners per unit length",
               title="Curvature density per component", legend=false)

    if outbase !== nothing
        savefig(plt1, outbase * "_hist_length.png")
        savefig(plt2, outbase * "_hist_tau.png")
        savefig(plt3, outbase * "_hist_curvdensity.png")
        savefig(plt4, outbase * "_scatter_len_tau.png")
        savefig(plt5, outbase * "_scatter_len_curvdensity.png")
        savefig(plt6, outbase * "_bar_comp_nstrands.png")
        savefig(plt7, outbase * "_bar_comp_totlen.png")
        savefig(plt8, outbase * "_bar_comp_meantau.png")
        savefig(plt9, outbase * "_bar_comp_curvdensity.png")
    end

    return (; plt1, plt2, plt3, plt4, plt5, plt6, plt7, plt8, plt9)
end

# === Run it now ===
strands_df, components_df = build_backbone_tables(strands, comps)

paths = save_backbone_metrics(strands_df, components_df;
    outdir = "logs/metrics",      # change if you want
    prefix = "backbone"           # <- your prefix goes here
)

plots = plot_backbone_metrics(strands_df, components_df; outbase = paths.base)
println("Saved tables to:\n  ", paths.strands_csv, "\n  ", paths.components_csv, "\n  ", paths.json)
println("Saved figures with base prefix:\n  ", paths.base, "_*.png")



# for i in 1:10
#     for j in i+1:10
#         if is_collision(x[i],y[1],x[j],y[j],1)
#             println("!")
#         end 
#     end
# end 

# generate_plots(state,config)