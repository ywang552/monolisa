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

generate_plots(state, config; output_prefix=prefix)

state.rotation

x = state.x_coords
y = state.y_coords;

# for i in 1:10
#     for j in i+1:10
#         if is_collision(x[i],y[1],x[j],y[j],1)
#             println("!")
#         end 
#     end
# end 

# generate_plots(state,config)