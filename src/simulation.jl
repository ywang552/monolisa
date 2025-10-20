# ==============================================
# File: simulation.jl
# Version: v1.0.0
# Description: Core simulation logic for MonoLisa, including monomer placement.
# Author: [Youchuan Wang]
# Date: [11/17/2024]
# ==============================================

include("dependencies.jl")
include("utils.jl")

const EPS = 1e-9  # put near top of file
const CONTACT_SCALE = 1.02  # keep a single source of truth
const TAU_ROW   = 0.05   # row threshold (knock out tiny mass)
const GAMMA_ROW = 2.0    # row sharpening (>=1)
const TOPK_ROWS = 40      # keep only top-k rows per parent (set 0 to disable)
const TAU_COL   = 0.03
const GAMMA_COL = 1.0


@inline dist2(x1,y1,x2,y2) = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)
# put near the top once
lo2(r::Real; scale::Real=CONTACT_SCALE) = (2r * (scale * 0.98))^2  # inner guard (optional)
hi2(r::Real; scale::Real=CONTACT_SCALE) = (2r * scale)^2           # must match find_contact_edges


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
    # last_to_check = 50

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

    # Reject if this box is full — do NOT relocate
    if boxNum[xtmp, ytmp] >= boxCap
        return false, -2   # "box full"
    end


    # # Ensure the box is within capacity
    # while boxNum[xtmp, ytmp] >= boxCap
    #     xtmp, ytmp = random_box_shift(xtmp, ytmp)  # Refactored into a helper function
    #     xtmp, ytmp = max(xtmp, 1), max(ytmp, 1)    # Ensure coordinates stay within bounds
    # end

    # # Adjust new_x and new_y based on box coordinates
    # new_x = (xtmp - 1) * boxSize + (new_x % boxSize)
    # new_y = (ytmp - 1) * boxSize + (new_y % boxSize)

    # Build column weights by multiplying rows
    # ensure we index with Ints (not Any)
    angles_for_pos = relative_angle[idx] isa AbstractVector ?
                    Int[relative_angle[idx]...] : [Int(relative_angle[idx])]

    col_weights = zeros(Float64, size(W, 2))
    for a in angles_for_pos
        v = @view(W[a, :])                           # <-- parenthesize the @view target
        col_weights .+= (max.(v .- TAU_COL, 0.0)).^(GAMMA_COL)
    end
    col_weights .= max.(col_weights, eps())
    ridx = sample(1:length(col_weights), Weights(col_weights))


    if all(iszero, col_weights)
        return false, -4   # invalid position: no joint probability
    end

    col_weights .= max.(col_weights, eps())
    ridx = sample(1:length(col_weights), Weights(col_weights))

    
    # 3×3 neighborhood sweep around candidate's box (xtmp, ytmp)
    sc = 0
    collision_detected = false
    new_NN_updates = Dict{Int,Int}()
    touch_contacts = Int[]               # raw list of all touched neighbors

    # clamp to grid bounds
    xlo = max(1, xtmp - 1); xhi = min(size(boxNum, 1), xtmp + 1)
    ylo = max(1, ytmp - 1); yhi = min(size(boxNum, 2), ytmp + 1)
    @inbounds for i1 in xlo:xhi
        for i2 in ylo:yhi
            inds = get(boxList, (i1, i2), Int[])
            for idx1 in inds
                d2 = dist2(new_x, new_y, x_coords[idx1], y_coords[idx1])

                if d2 < lo2(radius)                         # overlap => collision
                    collision_detected = true
                    break
                elseif d2 <= hi2(radius)                    # touching => valid neighbor
                    if NN[idx1] + 1 > NNrestriction
                        collision_detected = true
                        break
                    end
                    push!(touch_contacts, idx1)
                    new_NN_updates[idx1] = get(new_NN_updates, idx1, NN[idx1]) + 1
                    sc += 1
                end
            end
        end
        collision_detected && break
    end

    # (optional) cap the new monomer too
    if !collision_detected && sc > NNrestriction
        return false, -3
    end

    # Require at least one valid contact; otherwise reject
    # if !collision_detected && sc == 0
    #     return false, -6   # "no supporting neighbor"
    # end
    
    if collision_detected
        return false, -4
    else
        for (idx1, newval) in new_NN_updates
            NN[idx1] = newval
        end
    end


    # Add new monomer to the state
    new_index = size(x_coords, 1) + 1
    push!(x_coords, new_x)
    push!(y_coords, new_y)
    push!(NN, sc)

    # keep degree vector in sync
    if length(state.degree) < new_index
        resize!(state.degree, new_index)
    end
    state.degree[new_index] = 0


    # Update box data structures
    update_box_list!(boxList, boxNum, new_x, new_y, new_index, :add, config)

    # Update rotation and matrix K
    old_x, old_y = x_coords[tar[idx]], y_coords[tar[idx]]

    dx, dy = new_x - old_x, new_y - old_y
    tt_deg = round(Int, atan(-dy, -dx) * 180 / pi)  # cast to Int
    rot = mod(tt_deg - ridx * 5, 360)               # 0..359, stays Int
    push!(rotation, rot)

    for a in angles_for_pos
        K[a, ridx] += 1
    end

    unique!(touch_contacts)  # remove duplicates

    # always ensure parent is included
    parent = tar[idx]
    if parent != new_index && !(parent in touch_contacts)
        push!(touch_contacts, parent)
    end

    parent = tar[idx]
    @inbounds for v in touch_contacts
        v == new_index && continue
        u = v < new_index ? v : new_index
        w = v < new_index ? new_index : v
        push!(state.edges, (u, w))
        state.degree[v]        += 1
        state.degree[new_index] += 1
    end



    return true, boxIndx, [new_x, new_y]
end


# ============================
# Function: available_space(x_coords,
# Description: [Add a short description of what this function does.]
# Inputs: [List inputs here.]
# Outputs: [List outputs here.]
# ============================
# helper to delete multiple indices safely (reverse order)
function _delmany!(A::AbstractVector, idxs::AbstractVector{Int})
    for j in Iterators.reverse(sort(idxs))
        deleteat!(A, j)
    end
    return A
end

function _prune_parent_topk!(
    available_x::Vector{Float64},
    available_y::Vector{Float64},
    available_p::Vector{Float64},
    relative_angle::Vector{Any},
    target::Vector{Int},
    start_len::Int,             # length(available_p) BEFORE pushing this parent's angles
    K::Int
)
    K <= 0 && return
    lo = start_len + 1
    hi = length(available_p)
    lo > hi && return
    idxs = lo:hi
    length(idxs) <= K && return

    # find top-K (largest available_p) among this parent's pushes
    pool = collect(idxs)
    keep_perm = partialsortperm(pool, 1:K, by = j -> available_p[j], rev = true)
    keep = Set(pool[keep_perm])
    drop = [j for j in pool if j ∉ keep]
    isempty(drop) && return

    _delmany!(available_x, drop)
    _delmany!(available_y, drop)
    _delmany!(available_p, drop)
    _delmany!(relative_angle, drop)
    _delmany!(target, drop)
end

function available_space(x_coords, y_coords, radius, W, rotation, NN, LAST_TO_CHECK, NNrestriction)
    available_x = Float64[]           # candidate x
    available_y = Float64[]           # candidate y
    available_p = Float64[]           # per-candidate mass (pre-merge)
    relative_angle = Vector{Any}()    # row index or Vector{Int} of rows
    target = Int[]                    # parent index (proposing monomer)

    if length(x_coords) < LAST_TO_CHECK
        istart, iend = 1, length(x_coords)
    else
        istart, iend = length(x_coords) - LAST_TO_CHECK + 1, length(x_coords)
    end

    for i in istart:iend
        if NN[i] < NNrestriction
            avail = identify_available_angles(x_coords, y_coords, radius, i, NN, NNrestriction)

            # remember where this parent's block starts
            start_len = length(available_p)

            for angle in avail
                angleR   = calculate_relativeAngle(rotation[i], angle)
                angleidx = aidx(angleR)

                # sharpened row mass
                v = @view(W[angleidx, :])  # parenthesize target for @view
                p_raw   = sum(v)
                p_sharp = max(p_raw - TAU_ROW, 0.0)^GAMMA_ROW
                if p_sharp > 0
                    new_x = x_coords[i] + 2 * radius * cosd(angle)
                    new_y = y_coords[i] + 2 * radius * sind(angle)
                    push!(available_x, new_x)
                    push!(available_y, new_y)
                    push!(available_p, p_sharp)
                    push!(relative_angle, angleidx)
                    push!(target, i)
                end
            end

            # keep only top-K angles for this parent
            _prune_parent_topk!(available_x, available_y, available_p, relative_angle, target,
                                start_len, TOPK_ROWS)
        end
    end

    unique_x, unique_y, updated_probabilities, unique_angles, unique_targets =
        update_probabilities_fast(W, available_x, available_y, available_p, relative_angle, target)

    return unique_x, unique_y, updated_probabilities, unique_angles, unique_targets
end


"""
Group overlapping candidate positions and combine their probabilities.

- Overlap is decided by quantizing (x,y) to a grid set by `atol`.
- For a bucket with multiple candidates, we COMBINE by PRODUCT:
    p_pos = ∏_k p_k,  where  p_k = sum(W[row_k, :]).
- Supports the case where `relative_angle[i]` is a single Int or a Vector{Int}.
- Ensures strictly-positive outputs for Weights() by clamping with eps().

Returns:
    unique_x::Vector{Float64},
    unique_y::Vector{Float64},
    updated_probabilities::Vector{Float64},
    unique_angles::Vector{Vector{Int}},    # all rows (5° bins) that support each (x,y)
    unique_targets::Vector{Int}            # representative parent index for each (x,y)
"""
function update_probabilities_fast(
    W::AbstractMatrix,
    available_x::AbstractVector,
    available_y::AbstractVector,
    available_p::AbstractVector,
    relative_angle,     # Vector{Int} or Vector{Vector{Int}}
    target::AbstractVector;
    atol::Float64 = 0.05,    # positional merge tolerance (world units)
    rule::Symbol = :sharp_sum  # :sharp_sum | :sum | :product
)
    @assert length(available_x) == length(available_y) == length(available_p) ==
            length(relative_angle) == length(target)

    # Quantize (x,y) to merge near-identical positions (fast hash key)
    inv_step = 1 / atol
    key(x, y) = (round(x * inv_step) / inv_step, round(y * inv_step) / inv_step)

    # Bucket indices by position key
    groups = Dict{Tuple{Float64,Float64}, Vector{Int}}()
    for i in eachindex(available_x)
        k = key(available_x[i], available_y[i])
        push!(get!(groups, k, Int[]), i)
    end

    n = length(groups)
    unique_x              = Vector{Float64}(undef, n)
    unique_y              = Vector{Float64}(undef, n)
    unique_angles         = Vector{Vector{Int}}(undef, n)
    unique_targets        = Vector{Int}(undef, n)
    updated_probabilities = Vector{Float64}(undef, n)

    # Helper: combine probabilities in a bucket according to `rule`
    # p_k is based on row sums of W for that candidate's angle row(s).
    function calc_prob!(idxs::Vector{Int})::Float64
        if rule === :sum && length(idxs) > 1
            # legacy/fallback behavior
            return max(sum(@view available_p[idxs]), eps())
        end

        if rule === :sharp_sum
            acc = 0.0
            @inbounds for j in idxs
                rows_j = relative_angle[j]
                if rows_j isa AbstractVector
                    pj = 0.0
                    for r in rows_j
                        pj += max(sum(@view W[r, :]) - TAU_ROW, 0.0)^GAMMA_ROW
                    end
                    acc += pj
                else
                    r = rows_j::Int
                    acc += max(sum(@view W[r, :]) - TAU_ROW, 0.0)^GAMMA_ROW
                end
            end
            return max(acc^(1/GAMMA_ROW), eps())   # generalized mean (~softmax)
        end

        # PRODUCT of per-edge masses. A candidate may carry one row or many (Vector{Int})
        # For a candidate with multiple rows, we multiply their row-sums first.
        acc = 1.0
        for j in idxs
            rows_j = relative_angle[j]
            if rows_j isa AbstractVector
                pj = 1.0
                @inbounds for r in rows_j
                    pj *= max(sum(@view W[r, :]), eps())
                end
                acc *= pj
            else
                r = rows_j::Int
                acc *= max(sum(@view W[r, :]), eps())
            end
        end
        return max(acc, eps())
    end

    i = 1
    for (pos, idxs) in groups
        unique_x[i] = pos[1]
        unique_y[i] = pos[2]

        # collect/union all supporting angle rows for this position
        rows = Int[]
        for j in idxs
            rj = relative_angle[j]
            if rj isa AbstractVector
                append!(rows, rj)
            else
                push!(rows, rj)
            end
        end
        unique_angles[i] = unique!(rows)

        # pick a representative parent for bookkeeping (first is fine)
        unique_targets[i] = target[idxs[1]]

        # probability combine
        if length(idxs) == 1
            updated_probabilities[i] = max(available_p[idxs[1]], eps())
        else
            updated_probabilities[i] = calc_prob!(idxs)
        end
        i += 1
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
    seed_new_start_far!(state, config; dmin=100.0, dmax=500.0, tries=500)

Attempts to place a new monomer far from the existing structure.
- Chooses a random side of the current bounding box.
- Places it at a random distance between `dmin` and `dmax`.
- Checks grid capacity and collisions (using the same rules as `add_monomer`).
- If successful, commits the monomer with **no edges**.

Returns:
- `true` if a new seed was placed.
- `false` if all attempts failed.
"""
function seed_new_start_far!(state, config; dmin::Float64=100.0, dmax::Float64=500.0, tries::Int=500)
    xs, ys = state.x_coords, state.y_coords
    isempty(xs) && return false

    # Bounding box of current monomers
    xmin, xmax = minimum(xs), maximum(xs)
    ymin, ymax = minimum(ys), maximum(ys)

    radius  = state.radius
    boxNum  = state.boxNum
    boxList = state.boxList
    nrows, ncols = size(boxNum)

    @inbounds for _ in 1:tries
        # Random side and distance
        d    = rand()*(dmax - dmin) + dmin
        side = rand(1:4)

        x = 0.0; y = 0.0
        if side == 1
            x = xmin - d; y = rand()*(ymax - ymin) + ymin
        elseif side == 2
            x = xmax + d; y = rand()*(ymax - ymin) + ymin
        elseif side == 3
            y = ymin - d; x = rand()*(xmax - xmin) + xmin
        else
            y = ymax + d; x = rand()*(xmax - xmin) + xmin
        end

        # Map to grid box
        row, col = get_box_coordinates(x, y, config.box_size)
        if !(1 <= row <= nrows && 1 <= col <= ncols)
            continue
        end
        # Box capacity check
        if boxNum[row, col] >= config.box_capacity
            continue
        end

        # Collision check in 3×3 neighborhood (like add_monomer)
        collision = false
        xlo = max(1, row - 1); xhi = min(nrows, row + 1)
        ylo = max(1, col - 1); yhi = min(ncols, col + 1)
        for r in xlo:xhi
            for c in ylo:yhi
                inds = get(boxList, (r, c), Int[])
                for idx1 in inds
                    if dist2(x, y, xs[idx1], ys[idx1]) < lo2(radius)
                        collision = true
                        break
                    end
                end
                collision && break
            end
            collision && break
        end
        collision && continue

        # --- Commit seed (no edges) ---
        new_index = length(xs) + 1
        push!(state.x_coords, x)
        push!(state.y_coords, y)
        push!(state.NN, 0)
        if length(state.degree) < new_index
            resize!(state.degree, new_index)
        end
        state.degree[new_index] = 0
        push!(state.rotation, 0)  # can also randomize if desired

        update_box_list!(state.boxList, state.boxNum, x, y, new_index, :add, config)
        return true
    end

    return false
end


"""
Seed a short strand far from the existing structure and grow it in-place.

Returns: number of monomers placed in this seeding attempt (>=0).
"""
function seed_new_strand_far!(
    state, config;
    L::Int = state.last_to_check,    # target length of the new strand
    dmin::Float64 = 100.0,
    dmax::Float64 = 500.0,
    tries::Int = 500,                # spot attempts for the first seed
    per_node_fail::Int = 2000        # max add_monomer fails per subsequent node
)
    # 1) place the first seed far away
    ok = seed_new_start_far!(state, config; dmin=dmin, dmax=dmax, tries=tries)
    ok || return 0

    # indices bookkeeping
    first_idx = length(state.x_coords)          # the seed we just added
    placed = 1

    # 2) grow up to L-1 more nodes using the normal add_monomer pipeline,
    #    but anchor ONLY on the most recent nodes from this new strand.
    #    We achieve this by keeping last_to_check equal to how many nodes
    #    this new strand currently has.
    local_window = 1
    for _ in 2:L
        local_fails = 0
        # anchor only to the new strand's recent nodes
        state.last_to_check = local_window      # available_space uses this window
        # try until we place one, or give up for this node
        while true
            rs = add_monomer(state, config, 1)  # boix arg is unused upstream anyway
            if length(rs) == 3
                x = rs[1]
            else
                x = rs[1]
            end
            if x
                placed += 1
                local_window = min(local_window + 1, L)
                break
            else
                local_fails += 1
                if local_fails >= per_node_fail
                    # couldn't place this node; stop extending the seed
                    return placed
                end
            end
        end
    end
    return placed
end


"""
Runs the simulation and saves the minimal state.

Arguments:
- `fn`: Simulation identifier (if needed).
- `save_path`: File path to save the minimal state (default: "saved_states/minimal_state.jld2").

Returns:
- The final SimulationState.
"""
# function F(config; log_path=pwd()*"/"*"logs/", save_path=pwd()*"/"*"saved_states/minimal_states/")
#     # Initialize simulation state
#     state = initialize_simulation(config)
#     cc = 0
#     boix = 1
#     num_monomers = state.MN

#     log_path = log_path*"$(basename(config.file_path)[1:end-4]).log"
#     log_message(log_path, "Simulation started at: $(Dates.now())")
#     log_message(log_path, "Configuration: $(config)")

#     # Ensure the save folder exists
#     log_dir = dirname(log_path)
#     save_dir = dirname(save_path)
#     if !isdir(save_dir)
#         mkdir(save_dir)
#         println("Created save directory: $save_dir")
#         println("Created log directory: $log_dir")
#     end

#     # Set up progress bar
#     p = Progress(num_monomers, desc="Simulating Monomers", dt = config.prog)

#     try
#         while length(state.x_coords) < num_monomers
#             rs = add_monomer(state, config, boix)

#             if length(rs) == 3
#                 x, boix, npos = rs[1], rs[2], rs[3]
#             else
#                 x, boix = rs[1], rs[2]
#                 npos = -1
#             end
    
#             if x
#                 cc = 0
#                 ProgressMeter.next!(p)
#             else
#                 if boix == -1
#                     state.error_dic[1] += 1
#                     println("Break at ", size(state.x_coords, 1))
#                     break
#                 end
#                 if boix == -2
#                     cc += 1
#                     state.error_dic[2] += 1
#                 elseif boix == -3
#                     state.error_dic[3] += 1
#                 else
#                     cc += 1
#                     state.error_dic[4] += 1
#                 end
#             end
    
#             if x && (cc + 1) % 50 == 0
#                 println(size(state.x_coords, 1))
#                 println(state.error_dic)
#             end

    
#             if cc > 10000
#                 println("strand stalled (cc=$(cc)); seeding a new far mini-strand…")
#                 placed_seed = seed_new_strand_far!(
#                     state, config;
#                     L = state.last_to_check,      # e.g., 15 from your snippet
#                     dmin = 2., dmax = 4.,
#                     tries = 500, per_node_fail = 2000
#                 )
#                 if placed_seed > 0
#                     # keep anchoring within THIS new strand first
#                     state.last_to_check = min(state.last_to_check, placed_seed)
#                     cc = 0
#                     continue   # resume the SAME loop; normal growth takes over
#                 else
#                     println("no room to seed new mini-strand; breaking.")
#                     break
#                 end
#             end

#             if x && size(state.x_coords, 1) % 5000 == 0
#                 println("Placed ", size(state.x_coords, 1), " monomers")
#             end
#         end
#     catch e
#         threshold = 10  # Example: Minimum 3 monomers per grid region
#         # state = remove_small_islands!(state, threshold)
#         placed_monomer_number = length(state.x_coords)
#         pn = "$(placed_monomer_number)_$(basename(config.file_path)[1:end-4])"
#         # Save the current state as MinimalState on error
#         println("Error occurred: ", e)
#         log_message(log_path, "ERROR at $(Dates.now()): $e")
#         minimal_state = create_minimal_state(state)
#         save_minimal_state(save_path*"$(pn).jld2", minimal_state)
#         println("Minimal state saved after error.")
#         rethrow(e)  # Re-throw the error after saving
#     finally
#         threshold = 5  # Example: Minimum 3 monomers per grid region
#         # state = remove_small_islands!(state, threshold)
#         placed_monomer_number = length(state.x_coords)
#         pn = "$(placed_monomer_number)_$(basename(config.file_path)[1:end-4])"
#         # Save final results as MinimalState
#         minimal_state = create_minimal_state(state)
#         log_message(log_path, "Simulation completed at: $(Dates.now())")
#         save_minimal_state(save_path*"$(pn).jld2", minimal_state)
#         println("Final minimal state saved.")
#         # ProgressMeter.finish!(p)
#     end

#     println("Simulation complete!")
#     return state
# end

function F(config; log_path = joinpath(pwd(), "logs"),
                 save_path = joinpath(pwd(), "saved_states", "minimal_states"))

    # ===== Safenet knobs (tune as you like) =====
    target             = config.max_monomers
    stall_soft         = 5_000           # try local reseed after this many consecutive fails
    stall_hard         = 20_000          # escalate to far reseed/backoff
    max_global_plateau = 200_000         # last-ditch widen loop before admitting defeat
    reseed_tries       = 4               # how many far reseed attempts per escalation
    backoff_base       = (dmin=2.0, dmax=4.0)
    backoff_growth     = 1.5             # grow window on each escalation
    fatal_err_cap      = 5               # consecutive exceptions allowed

    # checkpoint cadence
    N_save       = 10_000                # save every N placements
    T_save_sec   = 300.0                 # or every 5 minutes
    last_save_ts = time()

    # ===== Init =====
    mkpath(save_path)
    mkpath(log_path)

    # derive filenames
    base = basename(config.file_path)
    stem = endswith(base, ".txt") ? base[1:end-4] : base
    log_file = joinpath(log_path, "$(stem).log")

    log_message(log_file, "Simulation started at: $(Dates.now())")
    log_message(log_file, "Configuration: $(config)")

    state = initialize_simulation(config)
    num_monomers = target

    # progress
    p = Progress(num_monomers, desc="Simulating Monomers", dt=config.prog)

    # internal counters
    stall_cc = 0                 # consecutive failures (no successful placement)
    hard_plateau = 0             # total failures since last big recovery
    fatal_errs = 0
    boix = 1

    # backoff window
    cur_dmin, cur_dmax = backoff_base

    # small helper: safe save checkpoint
    function save_checkpoint!(state; tag="chkpt")
        placed = length(state.x_coords)
        pn = "$(placed)_$(stem)_$(tag)"
        path = joinpath(save_path, pn * ".bin")

        open(path, "w") do io
            serialize(io, state)              # <-- FULL state
        end

        log_message(log_file,
            "Checkpoint saved at $(Dates.now()) (placed=$(placed), tag=$(tag), file=$(path))")
    end


    # small helper: staged reseed attempts with backoff
    function staged_reseed!(state, config; tries=reseed_tries, dmin=cur_dmin, dmax=cur_dmax)
        placed_seed = 0
        for t in 1:tries
            placed_seed = seed_new_strand_far!(
                state, config;
                L = state.last_to_check,
                dmin = dmin, dmax = dmax,
                tries = 1_000,           # be generous on search
                per_node_fail = 5_000
            )
            if placed_seed > 0
                state.last_to_check = min(state.last_to_check, placed_seed)
                return true
            end
        end
        return false
    end

    # ===== Main growth loop =====
    try
        while length(state.x_coords) < num_monomers
            # periodic autosave (by count or time)
            if (length(state.x_coords) % N_save == 0 && length(state.x_coords) > 0) ||
               (time() - last_save_ts >= T_save_sec)
                save_checkpoint!(state; tag="auto")
                last_save_ts = time()
            end

            # try to add one monomer
            rs = add_monomer(state, config, boix)

            # unpack result signature
            x::Bool = rs[1]
            boix    = rs[2]
            npos    = (length(rs) == 3) ? rs[3] : -1

            if x
                stall_cc = 0
                hard_plateau = 0
                ProgressMeter.next!(p)
                if length(state.x_coords) % 5_000 == 0
                    println("Placed ", length(state.x_coords), " monomers")
                    println("Errors: ", state.error_dic)
                end
                continue
            end

            # failure paths
            if boix == -1
                # previously: break (no room)
                # now: try staged recovery before giving up
                state.error_dic[1] += 1
                log_message(log_file, "boix=-1 encountered at $(length(state.x_coords)); attempting staged reseed…")

                # 1) local reseed with current window
                ok = staged_reseed!(state, config; dmin=cur_dmin, dmax=cur_dmax)
                if ok; stall_cc = 0; continue; end

                # 2) widen window and try again
                cur_dmin *= backoff_growth
                cur_dmax *= backoff_growth
                ok = staged_reseed!(state, config; dmin=cur_dmin, dmax=cur_dmax)
                if ok; stall_cc = 0; continue; end

                # 3) last-ditch: big window
                ok = staged_reseed!(state, config; dmin=cur_dmin*2, dmax=cur_dmax*2)
                if ok; stall_cc = 0; continue; end

                # If all fail, consider this a hard plateau step and keep looping until we exceed cap
                hard_plateau += 1_000
                if hard_plateau >= max_global_plateau
                    log_message(log_file, "Global plateau exceeded; saving and exiting loop.")
                    save_checkpoint!(state; tag="plateau")
                    break
                end

            elseif boix == -2
                stall_cc += 1
                state.error_dic[2] += 1

            elseif boix == -3
                state.error_dic[3] += 1
                stall_cc += 1
            else
                state.error_dic[4] += 1
                stall_cc += 1
            end

            # stall watchdogs
            if stall_cc >= stall_soft
                # Try a local reseed (no backoff change)
                ok = staged_reseed!(state, config; dmin=cur_dmin, dmax=cur_dmax)
                if ok
                    log_message(log_file, "Recovered from soft stall via local reseed.")
                    stall_cc = 0
                    continue
                end
            end

            if stall_cc >= stall_hard
                # escalate: widen window and try multiple reseeds
                cur_dmin *= backoff_growth
                cur_dmax *= backoff_growth
                ok = staged_reseed!(state, config; dmin=cur_dmin, dmax=cur_dmax)
                if ok
                    log_message(log_file, "Recovered from hard stall via widened reseed (dmin=$(cur_dmin), dmax=$(cur_dmax)).")
                    stall_cc = 0
                    continue
                end

                # if still stuck, record plateau (keeps loop alive but lets us exit eventually)
                hard_plateau += stall_cc
                stall_cc = 0
                if hard_plateau >= max_global_plateau
                    log_message(log_file, "Global plateau exceeded after hard stalls; saving and exiting loop.")
                    save_checkpoint!(state; tag="plateau")
                    break
                end
            end
        end

    catch e
        fatal_errs += 1
        println("Error occurred: ", e)
        log_message(log_file, "ERROR at $(Dates.now()): $e")
        save_checkpoint!(state; tag="error$(fatal_errs)")
        # do NOT rethrow; try to continue unless too many fatal errors
        if fatal_errs < fatal_err_cap
            @warn "Recovering from exception (attempt $fatal_errs of $fatal_err_cap) and continuing…"
            retry  # implicit by falling through to while; nothing special needed
        else
            @error "Too many fatal errors; exiting loop."
        end

    finally
        placed = length(state.x_coords)
        pn = "$(placed)_$(stem)"
        minimal_state = create_minimal_state(state)
        log_message(log_file, "Simulation completed at: $(Dates.now()) (placed=$(placed))")
        save_checkpoint!(state, tag = "final")
        try
            ProgressMeter.finish!(p)
        catch
            # ignore if already finished
        end
    end

    println("Simulation complete! placed=$(length(state.x_coords)) / target=$(target)")
    return state
end
