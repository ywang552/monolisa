using LinearAlgebra
using StatsBase
using NearestNeighbors
using VoronoiDelaunay
using Plots
using CSV
using DataFrames

# --- Load flat (x y x y x y ...) coordinate file ---
function load_flat_position_file(path::String)
    raw_line = read(path, String)
    nums = parse.(Float64, split(raw_line))
    @assert iseven(length(nums)) "Expected even number of values (x, y pairs)"
    coords = reshape(nums, 2, :)'  # N × 2
    return Matrix(Matrix(coords))
end

# --- Inter-monomer Distances ---
function pairwise_distances(positions::Matrix{Float64})
    tree = KDTree(positions')
    dists = []
    for i in 1:size(positions,1)
        idxs, ds = knn(tree, positions[i,:], 2)  # include self
        push!(dists, ds[2])  # second is nearest neighbor
    end
    return dists
end

# --- Voronoi Region Areas ---
function voronoi_area_distribution(positions::Matrix{Float64})
    vd = VD()
    for i in 1:size(positions, 1)
        push!(vd, positions[i,1], positions[i,2])
    end
    regions = voronoi(vd)
    areas = [area(r) for r in regions if isfinite(area(r))]
    return areas
end

# --- Grid Density ---
function local_density(positions::Matrix{Float64}, bin_size::Float64)
    xbins = floor(minimum(positions[:,1])):bin_size:ceil(maximum(positions[:,1]))
    ybins = floor(minimum(positions[:,2])):bin_size:ceil(maximum(positions[:,2]))
    counts = zeros(length(xbins)-1, length(ybins)-1)
    for i in 1:size(positions, 1)
        x, y = positions[i,1], positions[i,2]
        xi = searchsortedfirst(xbins, x) - 1
        yi = searchsortedfirst(ybins, y) - 1
        if xi ≥ 1 && yi ≥ 1 && xi ≤ size(counts,1) && yi ≤ size(counts,2)
            counts[xi, yi] += 1
        end
    end
    return counts
end

# --- Histogram Plot ---
function plot_histogram(data, title_str, nbins=50)
    histogram(data, bins=nbins, title=title_str, xlabel="Value", ylabel="Count")
end

function fast_pairwise_distances(coords::Matrix{Float64})
    # Input: coords is N×2, where each row is a monomer's (x, y)
    tree = KDTree(coords')                     # KDTree needs points as columns
    _, dists = knn(tree, coords', 2)           # Query ALL at once
    return [x[1] for x in dists]              # Take the 2nd closest (non-self)
end



function local_density(positions::Matrix{Float64}, bin_size::Float64)
    x_min, x_max = extrema(positions[:, 1])
    y_min, y_max = extrema(positions[:, 2])
    
    x_edges = x_min:bin_size:x_max
    y_edges = y_min:bin_size:y_max
    
    counts = zeros(Int, length(x_edges)-1, length(y_edges)-1)

    for i in 1:size(positions, 1)
        x, y = positions[i, :]
        xi = searchsortedlast(x_edges, x)
        yi = searchsortedlast(y_edges, y)
        if xi < size(counts, 1) + 1 && yi < size(counts, 2) + 1 && xi ≥ 1 && yi ≥ 1
            counts[xi, yi] += 1
        end
    end

    return counts, x_edges, y_edges
end


positions = load_flat_position_file("saved_states//minimal_states//hc2AFQ81NN139Q_Cpop_grids.txt")
x_min, x_max = extrema(positions[:, 1])
y_min, y_max = extrema(positions[:, 2])
println("X range: ", x_min, " to ", x_max)
println("Y range: ", y_min, " to ", y_max)


counts, xedges, yedges = local_density(positions, 10.0)
using Statistics

function print_density_stats(counts::Matrix{Int})
    flat_counts = vec(counts)  # Flatten 2D matrix into 1D vector of length 400

    # Only keep bins that have at least 1 monomer
    nz_counts = flat_counts[flat_counts .> 0]

    # --- Explanation for each statistic ---

    # Total number of bins = total spatial regions (e.g., 20×20 = 400)
    println("Total bins:         ", length(flat_counts))

    # Number of bins that have at least 1 monomer in them
    println("Non-zero bins:      ", count(>(0), flat_counts))  # e.g. 128 bins occupied

    # Sparsity: How much of the space is actually used
    # = (# of nonzero bins) / (total bins), as a percentage
    # Example: if 128 out of 400 bins are used => 32.0%
    println("Sparsity (%):       ", round(100 * count(>(0), flat_counts) / length(flat_counts), digits=2))

    # Mean number of monomers per bin, *including empty bins*
    # This gives an idea of how "dense" the simulation is overall
    println("Mean density:       ", round(mean(flat_counts), digits=2))

    # Standard deviation across all bins, including zeros
    # Indicates how unevenly the monomers are spread
    # A very high std (e.g. 2359.64) suggests some bins are extremely dense
    println("Std. dev (all bins):", round(std(flat_counts), digits=2))

    # Mean of only non-empty bins (ignores zeros)
    # Tells us how crowded occupied bins are
    # Example: if one bin has 46662 and the others ~10, this will be skewed high
    println("Mean (non-zero):    ", round(mean(nz_counts), digits=2))

    # Std dev of only non-empty bins
    # Gives variation among active regions only — a better idea of spread if we ignore empty space
    println("Std. dev (non-zero):", round(std(nz_counts), digits=2))

    # Bin with the most monomers — the densest cluster
    # This can signal a placement bug if it's way too high
    println("Max bin count:      ", maximum(flat_counts))

    # Least populated non-empty bin — smallest meaningful region
    println("Min (non-zero bin): ", minimum(nz_counts))
end

print_density_stats(counts)