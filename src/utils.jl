# ==============================================
# File: utils.jl
# Version: v1.0.0
# Description: Utility functions for MonoLisa, including collision detection and probability calculations.
# Author: [Youchuan Wang]
# Date: [11/17/2024]
# ==============================================

include("dependencies.jl")

"""
Function: is_collision
# Conditional check
Description: Determines if two monomers collide, based on their positions and radius.
# Conditional check
This is used to check if two monomers' spaces overlap, potentially indicating a collision.
Inputs:
    - x1, y1: Coordinates of the first monomer.
    - x2, y2: Coordinates of the second monomer.
    - radius: The radius of the monomers.
Returns:
    - A boolean indicating whether a collision occurs (True/False).
"""

# ============================
# Function: is_collision(x1,
# Description: [Add a short description of what this function does.]
# Inputs: [List inputs here.]
# Outputs: [List outputs here.]
# ============================
function is_collision(x1, y1, x2, y2, radius)
    distance_squared = (x1 - x2)^2 + (y1 - y2)^2
    return distance_squared < (2 * radius)^2  # Compare squared distance with squared diameter
end

"""
Function: is_neighbor
# Conditional check
Description: Determines if two monomers are neighbors, based on their coordinates.
# Conditional check
This is typically used to check if two monomers are adjacent in a grid system.
Inputs:
    - x1, y1: Coordinates of the first monomer.
    - x2, y2: Coordinates of the second monomer.
    - threshold: Maximum allowable distance between the monomers to be considered neighbors.
Returns:
    - A boolean indicating whether the two monomers are neighbors (True/False).
"""

# ============================
# Function: is_neighbor(x1,
# Description: [Add a short description of what this function does.]
# Inputs: [List inputs here.]
# Outputs: [List outputs here.]
# ============================
function is_neighbor(x1, y1, x2, y2, threshold)
    distance_squared = (x1 - x2)^2 + (y1 - y2)^2
    return distance_squared <= threshold^2  # Check if the distance is within the threshold
end

# ============================
# Function: calculate_relativeAngle(rotation,
# Description: [Add a short description of what this function does.]
# Inputs: [List inputs here.]
# Outputs: [List outputs here.]
# ============================
function calculate_relativeAngle(rotation, r)
    (r - rotation + 360) % 360
end 

# ============================
# Function: aidx(angle)
# Description: [Add a short description of what this function does.]
# Inputs: [List inputs here.]
# Outputs: [List outputs here.]
# ============================
function aidx(angle)
    if angle == 0 || angle == 360
        angleIdx = 1
    else 
        angleIdx = angle / 5 + 1
    end 

    angleIdx = Int(floor(angleIdx))
end 


# ============================
# Function: calculate_probabilities(W,
# Description: [Add a short description of what this function does.]
# Inputs: [List inputs here.]
# Outputs: [List outputs here.]
# ============================
function calculate_probabilities(W, indices, relative_angle)
    angles = relative_angle[indices]  # Get the relevant angles from relative_angle array
    return prod(W[angles, :])  # Multiplying the probabilities for the extracted angles
end


"""
Function: get_box_coordinates
# Loop through the specified range
Description: Calculates the box indices for given coordinates (x, y) based on box size.
Inputs:
    - x, y: Coordinates of the point.
    - box_size: Size of each box.
Returns:
    - A tuple (row, col) indicating the box indices.
"""
# ============================
# Function: get_box_coordinates(x::Float64,
# Description: [Add a short description of what this function does.]
# Inputs: [List inputs here.]
# Outputs: [List outputs here.]
# ============================
function get_box_coordinates(x::Float64, y::Float64, box_size::Int)
    row = Int((x - (x % box_size)) / box_size) + 1
    col = Int((y - (y % box_size)) / box_size) + 1
    return (row, col)
end

# ============================
# Function: update_box_list!(boxList::Dict{Tuple{Int,
# Description: [Add a short description of what this function does.]
# Inputs: [List inputs here.]
# Outputs: [List outputs here.]
# ============================
function update_box_list!(
    boxList::Dict{Tuple{Int, Int}, Vector{Int}}, 
    boxNum::SparseMatrixCSC{Int, Int}, 
    x::Float64, y::Float64, 
    index::Int, action::Symbol, 
    config::Config
)
    (row, col) = get_box_coordinates(x, y, config.box_size)
    
    if action == :add
        if haskey(boxList, (row, col))
            push!(boxList[(row, col)], index)
        else
            boxList[(row, col)] = [index]
        end
        boxNum[row, col] += 1
    elseif action == :remove
        if haskey(boxList, (row, col))
            filter!(i -> i != index, boxList[(row, col)])
        end
        boxNum[row, col] -= 1
    else
        error("Invalid action. Use :add or :remove.")
    end
end


function random_box_shift(xtmp::Int, ytmp::Int)
    # Map random values to coordinate changes
    directions = Dict(
        1 => (1, 0),   # Right
        2 => (1, 1),   # Top-right
        3 => (0, 1),   # Up
        4 => (-1, 1),  # Top-left
        5 => (-1, 0),  # Left
        6 => (-1, -1), # Bottom-left
        7 => (0, -1),  # Down
        8 => (1, -1)   # Bottom-right
    )
    
    # Get the shift for the chosen direction
    value = rand(1:8)
    shift = directions[value]
    
    # Return updated coordinates
    return xtmp + shift[1], ytmp + shift[2]
end

"""
Saves a MinimalState to a file.

Arguments:
- `file_path`: The file path where the MinimalState will be saved.
- `minimal_state`: The MinimalState struct to save.
"""
function save_minimal_state(file_path::String, minimal_state::MinimalState)
    # Save the essential fields of MinimalState
    @save file_path minimal_state
    println("Minimal state saved to: $file_path")
end


"""
Removes small islands (grid regions with fewer monomers than a threshold).

Arguments:
- `state`: The SimulationState object containing the simulation data.
- `threshold`: Minimum number of monomers required in a grid region to retain it.

Returns:
- Updated SimulationState with small islands removed.
"""
function remove_small_islands!(state, threshold)
    removed_indices = Set{Int}()

    # Identify grid regions with fewer than the threshold
    for ((row, col), monomer_indices) in state.boxList
        if state.boxNum[row, col] < threshold
            union!(removed_indices, monomer_indices)  # Collect indices to remove
        end
    end

    # Convert the set of indices to a sorted vector
    removed_indices = sort(collect(removed_indices))

    # Build a mapping from old indices to new indices
    total_monomers = length(state.x_coords)
    remaining_indices = setdiff(1:total_monomers, removed_indices)
    index_map = Dict(old => new for (new, old) in enumerate(remaining_indices))

    # Batch remove monomers from the state
    deleteat!(state.x_coords, removed_indices)
    deleteat!(state.y_coords, removed_indices)
    deleteat!(state.NN, removed_indices)

    # Update boxNum and boxList
    new_boxList = Dict{Tuple{Int, Int}, Vector{Int}}()
    for ((row, col), monomer_indices) in state.boxList
        # Update indices to match the new index mapping
        updated_indices = filter(x -> haskey(index_map, x), monomer_indices)
        updated_indices = map(x -> index_map[x], updated_indices)

        if !isempty(updated_indices)
            new_boxList[(row, col)] = updated_indices
        else
            state.boxNum[row, col] = 0
        end
    end
    state.boxList = new_boxList

    println("Removed $(length(removed_indices)) monomers from small islands.")
    return state
end

"""
Logs a message to a log file.

Arguments:
- `log_file`: The file path to the log file.
- `message`: The message to log.
"""
function log_message(log_file::String, message::String)
    open(log_file, "a") do io
        println(io, message)
    end
end
