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
include("structure_analysis.jl")
include("longest_strand.jl")
include("faces.jl")


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
    box_capacity=10,
    monomer_radius=1,
    grid_size=2000,
    max_monomers=1000,
    # max_monomers=100,
    # file_path=ARGS[1],
    file_path="Claudins/wt2_newsep.txt",
    grid_overlay = false,
    prog = 1
)

function main()
    println("starting the placement...")
    # Pass configuration to simulation
    state = F(config, log_path = "large_strand\\logs\\", save_path = "large_strand\\placements\\")

    # Save results or generate plots
    return state
end


safe_stamp() = Dates.format(now(), "yyyymmdd-HHMMSS-sss")
stamp = safe_stamp()
state = main();





















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
