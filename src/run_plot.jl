using JLD2
using FilePathsBase  # For handling paths

# Define the path to the saved states folder and output folder
const SAVED_STATES_PATH = "saved_states/minimal_states/"
const OUTPUT_FOLDER = "plots/tmp"

# Ensure the output folder exists
isdir(OUTPUT_FOLDER) || mkdir(OUTPUT_FOLDER)

# Function to load all saved states from .jld2 files
function load_saved_states(path::String)
    @load path k 
    return k
end

fn = "89816_c4phoFelix_grids.jld2"
@load SAVED_STATES_PATH*fn minimal_state


config.file_path = split(config.file_path, "/")[1] * "/" * split(fn, "_")[2]

generate_plots(minimal_state, config, pwd()*"\\plots\\tmp\\")
