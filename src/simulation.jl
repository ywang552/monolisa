# ==============================================
# File: simulation.jl
# Version: v1.0.0
# Description: Core simulation logic for MonoLisa, including monomer placement.
# Author: [Youchuan Wang]
# Date: [11/17/2024]
# ==============================================

include("dependencies.jl")
include("utils.jl")

# ============================
# Function: add_monomer(x_coords,
# Description: [Add a short description of what this function does.]
# Inputs: [List inputs here.]
# Outputs: [List outputs here.]
# ============================
function add_monomer(state::SimulationState, config::Config, boxIndx::Int)
    # Extract parameters from state and config
    x_coords, y_coords, rotation, NN = state.x_coords, state.y_coords, state.rotation, state.NN
    radius, boxNum, boxList, K = state.radius, state.boxNum, state.boxList, state.K
    W = state.W
    last_to_check = state.last_to_check

    boxSize, NNrestriction, boxCap = config.box_size, config.nn_restriction, config.box_capacity

    # Get available positions
    available_x, available_y, available_p, relative_angle, tar =
        available_space(x_coords, y_coords, radius, W, rotation, NN, last_to_check, NNrestriction)

    if isempty(available_x)
        println("No available space to add a monomer.")
        return false, -1
    end

    # Choose a random position based on probabilities
    idx = sample(1:length(available_p), Weights(Float64.(available_p)))
    new_x, new_y = available_x[idx], available_y[idx]

    # Determine box coordinates for the new position
    (xtmp, ytmp) = get_box_coordinates(new_x, new_y, boxSize)

    # Ensure the box is within capacity
    while boxNum[xtmp, ytmp] >= boxCap
        xtmp, ytmp = random_box_shift(xtmp, ytmp)  # Refactored into a helper function
        xtmp, ytmp = max(xtmp, 1), max(ytmp, 1)    # Ensure coordinates stay within bounds
    end

    # Adjust new_x and new_y based on box coordinates
    new_x = (xtmp - 1) * boxSize + (new_x % boxSize)
    new_y = (ytmp - 1) * boxSize + (new_y % boxSize)

    # Update probabilities and angles
    angle = relative_angle[idx]
    w = reshape(sum(W[angle, :], dims=1), 72)
    ridx = sample(1:length(w), Weights(w))

    # Check for collisions and update neighbors
    sc = 0
    collision_detected = false
    new_NN_updates = Dict()

    # for i1 in max(1, xtmp - 1):xtmp + 1
    #     for i2 in max(1, ytmp - 1):ytmp + 1
    #         for i in 1:boxNum[i1, i2]
    #             idx1 = boxList[i1, i2][i]

    #             if is_neighbor(new_x, new_y, x_coords[idx1], y_coords[idx1], radius)
    #                 if is_collision(new_x, new_y, x_coords[idx1], y_coords[idx1], radius) || NN[idx1] + 1 > NNrestriction
    #                     collision_detected = true
    #                 else
    #                     new_NN_updates[idx1] = get(new_NN_updates, idx1, NN[idx1]) + 1
    #                     sc += 1
    #                 end
    #             end
    #         end
    #     end
    # end

    for i in eachindex(x_coords)
        if is_collision(new_x, new_y, x_coords[i], y_coords[i], radius)
            collision_detected = true
        end 
    end 

    if collision_detected
        return false, -4
    else
        for (idx, new_val) in new_NN_updates
            NN[idx] = new_val
        end
    end

    # Add new monomer to the state
    new_index = size(x_coords, 1) + 1
    push!(x_coords, new_x)
    push!(y_coords, new_y)
    push!(NN, sc)

    # Update box data structures
    update_box_list!(boxList, boxNum, new_x, new_y, new_index, :add, config)

    # Update rotation and matrix K
    old_x, old_y = x_coords[tar[idx]], y_coords[tar[idx]]
    dx, dy = new_x - old_x, new_y - old_y
    tt = round(atan(-dy, -dx) * 180 / pi)
    push!(rotation, (tt - ridx * 5 + 720) % 360)
    K[angle, ridx] .+= 1

    return true, boxIndx, [new_x, new_y]
end


# ============================
# Function: available_space(x_coords,
# Description: [Add a short description of what this function does.]
# Inputs: [List inputs here.]
# Outputs: [List outputs here.]
# ============================
function available_space(x_coords, y_coords, radius, W, rotation, NN, LAST_TO_CHECK, NNrestriction)
    available_x = []
    available_y = []
    available_p = []
    relative_angle = []
    target = []
    tmp = []
    if (length(x_coords) < LAST_TO_CHECK)
        for i in 1:length(x_coords)
            if(NN[i] < NNrestriction)
                available_angles = identify_available_angles(x_coords, y_coords, radius,i, NN, NNrestriction)
                for angle in 0:5:355
                    angleR = calculate_relativeAngle(rotation[i], angle)
                    angleidx1 = aidx(angleR)
                    p = sum(W[angleidx1, :])
                    if(p!=0)
                        new_x = round(x_coords[i] + 2 * radius * cosd(angle), digits = 3)
                        new_y = round(y_coords[i] + 2 * radius * sind(angle), digits = 3)
                        push!(available_x, new_x)
                        push!(available_y, new_y)
                        push!(available_p, p)
                        push!(relative_angle, angleidx1)
                        push!(target, i)
                    end 

                end
            end 
        end
    else
        for i in 1+length(x_coords)-LAST_TO_CHECK:length(x_coords)
            if(NN[i] < NNrestriction)
                available_angles = identify_available_angles(x_coords, y_coords, radius,i, NN, NNrestriction)
                for angle in available_angles
                    angleR = calculate_relativeAngle(rotation[i], angle)
                    angleidx1 = aidx(angleR)
                    p = sum(W[angleidx1, :])

                    if(p!=0)
                        new_x = round(x_coords[i] + 2 * radius * cosd(angle), digits = 3)
                        new_y = round(y_coords[i] + 2 * radius * sind(angle), digits = 3)
                        push!(available_x, new_x)
                        push!(available_y, new_y)
                        push!(available_p, p)
                        push!(relative_angle, angleidx1)
                        push!(target, i)
                    end 
            
                end
            end 
        end
    end

    overlapping_indices = []
    seen_positions = []  # To track seen (x, y) positions

    for i in 1:length(available_x)
        position = (available_x[i], available_y[i])
        found_overlap = false

        for (j, seen_position) in enumerate(seen_positions)
            if isapprox(seen_position[1], position[1]) && isapprox(seen_position[2], position[2])
                push!(overlapping_indices, i)  # Record the current index
                push!(overlapping_indices, j)    # Record the index of the previously seen position
                found_overlap = true
                break  # Exit inner loop since we found an overlap
            end
        end

        if !found_overlap
            push!(seen_positions, position)
        end
    end

    unique_x, unique_y, updated_probabilities, unique_angles, unique_targets = 
    update_probabilities_fast(W, available_x, available_y, available_p, relative_angle, target)

    unique_x, unique_y, updated_probabilities, unique_angles, unique_targets
end

# ============================
# Function: update_probabilities_fast(W,
# Description: [Add a short description of what this function does.]
# Inputs: [List inputs here.]
# Outputs: [List outputs here.]
# ============================
function update_probabilities_fast(W, available_x, available_y, available_p, relative_angle, target)
    overlapping_indices = Dict{Tuple{Float64, Float64}, Vector{Int}}()

    for i in 1:length(available_x)
        position = (available_x[i], available_y[i])

        found_overlap = false
        for (key, indices) in overlapping_indices
            if isapprox(key[1], position[1], atol=0.05) && isapprox(key[2], position[2], atol=0.05)
                push!(indices, i)
                found_overlap = true
                break
            end

        end

        if !found_overlap
            overlapping_indices[position] = [i]
        end
    end

    num_unique = length(overlapping_indices)
    unique_x = Vector{Float64}(undef, num_unique)
    unique_y = Vector{Float64}(undef, num_unique)
    unique_angles = Vector{Vector{Int}}(undef, num_unique)
    unique_targets = Vector{Int}(undef, num_unique)
    updated_probabilities = Vector{Float64}(undef, num_unique)

    index = 1
    for (position, indices) in overlapping_indices
        unique_x[index] = position[1]
        unique_y[index] = position[2]

        if length(indices) > 1
            new_prob = calculate_probabilities(W, indices, relative_angle)
            combined_angles = relative_angle[indices]  # Combine the angles for this group
            unique_angles[index] = combined_angles
            updated_probabilities[index] = new_prob
        else
            idx = indices[1]
            unique_angles[index] = [relative_angle[idx]]
            updated_probabilities[index] = available_p[idx]
        end

        unique_targets[index] = target[indices[1]]
        index += 1
    end

    return unique_x, unique_y, updated_probabilities, unique_angles, unique_targets
end

# ============================
# Function: identify_available_angles(x_coords,
# Description: [Add a short description of what this function does.]
# Inputs: [List inputs here.]
# Outputs: [List outputs here.]
# ============================
function identify_available_angles(x_coords, y_coords, radius, idx, NN, NNrestriction)
    angles = 0:5:355  # Angles at 5-degree intervals
    available_angles = Vector{Float64}()
    
    for angle in angles
        col = false

    
        if(!col)
            push!(available_angles,angle)            
        end

    end 

    return available_angles
end


"""
Runs the simulation and saves the minimal state.

Arguments:
- `fn`: Simulation identifier (if needed).
- `save_path`: File path to save the minimal state (default: "saved_states/minimal_state.jld2").

Returns:
- The final SimulationState.
"""
function run(fn; log_path=pwd()*"/"*"logs/", save_path=pwd()*"/"*"saved_states/minimal_states/")
    # Initialize simulation state
    state = initialize_simulation(config)

    cc = 0
    boix = 1
    num_monomers = state.MN

    log_path = log_path*"$(basename(config.file_path)[1:end-4]).log"
    log_message(log_path, "Simulation started at: $(Dates.now())")
    log_message(log_path, "Configuration: $(config)")

    # Ensure the save folder exists
    log_dir = dirname(log_path)
    save_dir = dirname(save_path)
    if !isdir(save_dir)
        mkdir(save_dir)
        println("Created save directory: $save_dir")
        println("Created log directory: $log_dir")
    end

    # Set up progress bar
    p = Progress(num_monomers, desc="Simulating Monomers", dt = config.prog)

    try
        while length(state.x_coords) < num_monomers
            rs = add_monomer(state, config, boix)

            if length(rs) == 3
                x, boix, npos = rs[1], rs[2], rs[3]
            else
                x, boix = rs[1], rs[2]
                npos = -1
            end
    
            if x
                cc = 0
                ProgressMeter.next!(p)
            else
                if boix == -1
                    state.error_dic[1] += 1
                    println("Break at ", size(state.x_coords, 1))
                    break
                end
                if boix == -2
                    cc += 1
                    state.error_dic[2] += 1
                elseif boix == -3
                    state.error_dic[3] += 1
                else
                    cc += 1
                    state.error_dic[4] += 1
                end
            end
    
            if x && (cc + 1) % 50 == 0
                println(size(state.x_coords, 1))
                println(state.error_dic)
            end
    
            if cc > 100
                println("broke: too many collisions/empty row")
                break
            end
    
            if x && size(state.x_coords, 1) % 5000 == 0
                println("Placed ", size(state.x_coords, 1), " monomers")
            end
        end
    catch e
        threshold = 10  # Example: Minimum 3 monomers per grid region
        state = remove_small_islands!(state, threshold)
        placed_monomer_number = length(state.x_coords)
        pn = "$(placed_monomer_number)_$(basename(config.file_path)[1:end-4])"
        # Save the current state as MinimalState on error
        println("Error occurred: ", e)
        log_message(log_path, "ERROR at $(Dates.now()): $e")
        minimal_state = create_minimal_state(state)
        save_minimal_state(save_path*"$(pn).jld2", minimal_state)
        println("Minimal state saved after error.")
        rethrow(e)  # Re-throw the error after saving
    finally
        threshold = 5  # Example: Minimum 3 monomers per grid region
        state = remove_small_islands!(state, threshold)
        placed_monomer_number = length(state.x_coords)
        pn = "$(placed_monomer_number)_$(basename(config.file_path)[1:end-4])"
        # Save final results as MinimalState
        minimal_state = create_minimal_state(state)
        log_message(log_path, "Simulation completed at: $(Dates.now())")
        save_minimal_state(save_path*"$(pn).jld2", minimal_state)
        println("Final minimal state saved.")
        # ProgressMeter.finish!(p)
    end

    println("Simulation complete!")
    return state
end