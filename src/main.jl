# ==============================================
# File: main.jl
# Version: v1.0.0
# Description: Entry point for MonoLisa simulation.
# Author: [Youchuan Wang]
# Date: [11/17/2024]
# ==============================================

include("dependencies.jl")
include("initialization.jl")
# include("simulation.jl")
include("simulation_forced_flag.jl")
include("utils.jl")
include("plot_results.jl")
include("structure_analysis.jl")
# include("longest_strand.jl")
include("compute_backbone_fast.jl")
include("faces.jl")
include("x.jl")
include("y.jl")
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
fn = "wt2_newsep"
config = Config(
    ft=0.0005,
    box_size=10,
    nn_restriction=3,
    box_capacity=15,
    monomer_radius=1,
    grid_size=2000,
    max_monomers=500,
    # max_monomers=100,
    # file_path=ARGS[1],
    file_path="Claudins/$(fn).txt",
    grid_overlay = false,
    prog = 1
)
function main(; s = "")
    println("starting the placement...")
    # Pass configuration to simulation
    state = F(config; log_path = pwd()*"\\large_strand\\logs\\", save_path = pwd()*"\\large_strand\\placements\\", stamp = s)

    # Save results or generate plots
    return state
end



"Return a short, filesystem-safe hash timestamp."
@inline function safe_stamp()
    t = string(Dates.now())
    h = bytes2hex(sha1(t))
    return h[1:10]     # 10-char hash prefix, plenty unique per second
end

stamp = safe_stamp()
state = main(;s = stamp);
# state = deserialize("large_strand\\placements\\8000_wt2_newsep_final.bin")

backbones = compute_backbone(state; λ=0.6, mode=:geodesic)
out_dir = joinpath(joinpath("plots", fn))
if !ispath(out_dir)
    mkdir(out_dir)
end 
prefix = joinpath(out_dir, fn)
generate_plots(state, config; bbs = backbones, output_prefix = prefix*"_$(stamp)", show_contour=true, tm_style=:nothing)






# Folder with your .bin states
data_dir = joinpath("large_strand","placements")
# Customize how filenames map to claudin labels, if needed
patterns = Dict(
    # "C4" => "C4",
    # "C2" => "C2",
    "hc5AF" => "C5",
    "wt2_newsep" => "C2",
    "hc15" => "C15",
    "c4_7" => "C4",

    # add more keys if filenames use other hints, e.g. "claudin4" => "C4"
)

out_dir = joinpath("plots", "cdf")
if !ispath(out_dir)
    mkdir(out_dir)
end 
area_floor  = 8.
groups = load_states_grouped(data_dir; patterns=patterns)

plots_hist = plot_hist_longest_by_type(groups; λ=0.6, mode=:geodesic, density=true)
ps = plot_histograms_by_claudin(data_dir, patterns, density = false, area_floor = area_floor)

for p in keys(ps)
    local prefix = joinpath(out_dir, "$(p)")
    savefig(ps[p], "$(prefix)_MeshAreaHistogram.png")
    savefig(plots_hist[p], "$(prefix)_LongestStrandHistogram.png")
end 

plt = plot_cdf_clusters_by_type(groups;
    colors=TRUE_COLORS,
    area_floor=8.2,
    use_ccdf=true,
    focus_top=(0.7, 1.0),
)
prefix = joinpath(out_dir, "mesh_area_cdf_zoomed.png")
savefig(plt, prefix)

plt = plot_cdf_clusters_by_type(groups;
    colors=TRUE_COLORS,
    area_floor=8.2,
    use_ccdf=true,
    focus_top=(0.0, 1.0), 
)
prefix = joinpath(out_dir, "mesh_area_cdf.png")
savefig(plt, prefix)
plt_cdf = plot_cdf_longest_by_type(groups; λ=0.6, mode=:geodesic, use_ccdf=false)

prefix = joinpath(out_dir, "longest_strand_cdf.png")
savefig(plt_cdf, prefix)


# """
#     report_overconnected(state; max_deg=3, use_geometric=false, contact_scale=1.02)

# Check each monomer's number of neighbors and report any with degree > max_deg.
# """
# function report_overconnected(state;
#     max_deg::Int=2,
#     use_geometric::Bool=false,
#     contact_scale::Real=1.02
# )

#     x, y, r = state.x_coords, state.y_coords, state.radius
#     xmin, xmax = minimum(x), maximum(x)
#     ymin, ymax = minimum(y), maximum(y)
#     pad = 2r
#     xlim = (xmin - pad, xmax + pad)
#     ylim = (ymin - pad, ymax + pad)
#     edge_list = find_contact_edges(x, y, r; xlim=xlim, ylim=ylim,
#                                     contact_scale=contact_scale)


#     # --- degree map ---
#     deg = Dict{Int,Int}()
#     for (u,v) in edge_list
#         u == v && continue
#         deg[u] = get(deg,u,0) + 1
#         deg[v] = get(deg,v,0) + 1
#     end

#     offenders = [i for (i,d) in deg if d > max_deg]

#     println("Neighbor-degree check (threshold = $max_deg)")
#     println("Total monomers: ", length(state.x_coords),
#             " | edges considered: ", length(edge_list))
#     if isempty(offenders)
#         println("✓ No monomers exceed degree $max_deg.")
#     else
#         println("⚠ ", length(offenders), " monomer(s) exceed degree $max_deg:")
#         for i in offenders
#             println("  monomer $i → degree $(deg[i])")
#         end
#     end

#     return (; max_deg, offenders, degrees=deg, edge_list)
# end


# ms = create_minimal_state(state)
# # # ms.edges
# # # state.edges
# # # state.last_to_check


# # # # print(length(state.x_coords))
# out_dir = joinpath(pwd(), "plots", "tmp")
# mkpath(out_dir)

# name_noext, _ = splitext(basename(config.file_path))
# prefix = joinpath(out_dir, name_noext*"_123")
# overlay_path = prefix * "overlay_curvature.png"
# # # println(typeof(state2))
# # # state2.edges
# # # state.edges
# backbones = compute_backbone(state; λ=0.6, mode=:geodesic)
# generate_plots(state, config; bbs = backbones, output_prefix = prefix*"arcs_$(stamp)", show_contour=true, tm_style=:nothing)
# # # state2

# x = state.x_coords
# y = state.y_coords;
# state.edges
# using CSV, DataFrames

# df_a = DataFrame(x = state.x_coords, y = state.y_coords, r = state.rotation)
# CSV.write("3D_info.csv", df_a)
