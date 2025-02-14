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
    grid_size=8000,
    # max_monomers=1000000,
    max_monomers=20000,
    # file_path=ARGS[1],
    file_path="Claudins/Cpop_grids 2.txt",
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

generate_plots(state, config, pwd()*"\\plots\\tmp\\"*basename(config.file_path)[1:end-4]*"_123")
