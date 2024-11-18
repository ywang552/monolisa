# ==============================================
# File: initialization.jl
# Version: v1.0.0
# Description: Handles simulation state initialization and data loading.
# Author: [Youchuan Wang]
# Date: [11/17/2024]
# ==============================================
include("dependencies.jl")


abstract type AbstractState end


mutable struct Config
    ft::Float64                # Threshold for probabilities
    box_size::Int              # Size of each grid box
    nn_restriction::Int        # Neighbor restriction
    box_capacity::Int          # Maximum monomers per box
    monomer_radius::Float64    # Radius of each monomer
    grid_size::Int             # Grid size
    max_monomers::Int          # Maximum number of monomers
    file_path::String          # Path to the input file
    grid_overlay::Bool
    # Constructor to allow keyword arguments
    function Config(; ft=0.0005, box_size=10, nn_restriction=3, box_capacity=16,
                     monomer_radius=1.0, grid_size=6000, max_monomers=100, file_path="8Cpop_grids.txt", grid_overlay = false)
        new(ft, box_size, nn_restriction, box_capacity, monomer_radius, grid_size, max_monomers, file_path, grid_overlay)
    end
end

config = Config()  # Uses default values




mutable struct SimulationState <: AbstractState
    x_coords::Vector{Float64}
    y_coords::Vector{Float64}
    rotation::Vector{Int}
    NN::Vector{Int}
    radius::Float64
    last_to_check::Int
    NNrestriction::Int
    boxNum::SparseMatrixCSC{Int, Int}
    boxList::Dict{Tuple{Int, Int}, Vector{Int}}
    boxCap::Int
    box_size::Int
    MN::Int
    W::Matrix{Float64}
    K::Matrix{Float64}
    error_dic::Vector{Int}
end

struct MinimalState <: AbstractState
    x_coords::Vector{Float64}
    y_coords::Vector{Float64}
    box_size::Int
    W::Matrix{Float64}
    K::Matrix{Float64}
    
end

function create_minimal_state(state::SimulationState)
    return MinimalState(
        state.x_coords,
        state.y_coords,
        state.box_size,
        state.W,
        state.K
    )
end


"""
Function: initialize_simulation
Description: Sets up all simulation variables and processes the W matrix.
Inputs:
    - MN: Number of monomers.
    - DS: Grid size.
    - boxCap: Maximum monomers per grid box.
    - boxSize: Size of each box in the grid.
    - NNrestriction: Max number of neighbors for a monomer.
    - radius: Radius of each monomer.
    - file_path: Path to the file containing the W matrix.
Outputs:
    - A dictionary containing all initialized variables.
"""

function initialize_simulation(config::Config)
    # Extract parameters from config
    MN = config.max_monomers
    DS = config.grid_size
    boxCap = config.box_capacity
    boxSize = config.box_size
    NNrestriction = config.nn_restriction
    radius = config.monomer_radius
    file_path = config.file_path

    # Initialize state
    x_coords = Float64[]
    y_coords = Float64[]
    rotation = [Int(18)]
    NN = [0]
    boxNum = spzeros(Int, DS, DS)
    boxList = Dict{Tuple{Int, Int}, Vector{Int}}()

    # Place the first monomer
    xsp = Int(DS * 5)
    ysp = Int(DS * 5)
    boxList[Int(xsp / boxSize) + 1, Int(ysp / boxSize) + 1] = [1]
    boxNum[Int(xsp / boxSize) + 1, Int(ysp / boxSize) + 1] = 1
    push!(x_coords, xsp)
    push!(y_coords, ysp)

    # Load and preprocess W matrix
    W = load_w_matrix(file_path)
    K = zeros(72, 72)

    # Return initialized state
    return SimulationState(
        x_coords, y_coords, rotation,
        NN, radius, 5,  # last_to_check
        NNrestriction, boxNum, boxList,
        boxCap, boxSize, MN, W, K,
        zeros(Int, 4)  # error_dic
    )
    
end

"""
Function: load_w_matrix
Description: Reads and processes the W matrix from a file.
Inputs:
    - file_path: Path to the file containing the W matrix data.
Outputs:
    - A normalized W matrix.
"""

"""
Function: load_w_matrix
Description: Reads and processes the W matrix from a file, ensuring proper dimensions.
Inputs:
    - file_path: Path to the file containing the W matrix data.
Outputs:
    - A normalized W matrix with dimensions 360x360.
"""
function load_w_matrix(file_path)
    # Read W matrix as a vector
    file = open(file_path, "r")
    raw_data = Float64[]  # Store all the values in a flat array
    for line in eachline(file)
        row = parse.(Float64, split(line))  # Parse each line into numbers
        append!(raw_data, row)
    end
    close(file)

    # Ensure the W matrix has the expected size
    if length(raw_data) != 360 * 360
        error("The W matrix file does not contain the expected number of elements (360x360).")
    end

    # Reshape into 360x360 matrix
    P_original = reshape(raw_data, 360, 360)

    # Downsample W matrix to 72x72 and normalize
    new_size = 72
    ratio = size(P_original, 1) ÷ new_size
    P_shrunk = zeros(new_size, new_size)

    for i in 1:new_size
        for j in 1:new_size
            start_row = 1 + (i - 1) * ratio
            end_row = i * ratio
            start_col = 1 + (j - 1) * ratio
            end_col = j * ratio

            # Extract the submatrix
            submatrix = P_original[start_row:end_row, start_col:end_col]
            P_shrunk[i, j] = mean(submatrix)
        end
    end

    # Normalize the downsampled matrix
    W = P_shrunk ./ sum(P_shrunk)
    msk2 = W .<= 0.0005
    W[msk2] .= 0
    return (W .- minimum(W)) ./ (maximum(W) - minimum(W))
end



