# === includes (ensure SimulationState & get_areas are defined before loading) ===
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

using Glob
# --- tiny loader for .bin files produced by Serialization.serialize ---
load_state(path::AbstractString) = open(path, "r") do io
    deserialize(io)  # expects a SimulationState defined by your includes above
end

# === config ===
data_dir = joinpath(pwd(), "large_strand", "placements")
out_dir  = joinpath(pwd(), "plots", "cdf")

patterns = Dict(
    "hc5AF"      => "C5",
    "wt2_newsep" => "C2",
    "hc15"       => "C15",
    "c4_7"       => "C4",
)

max_files = 5

# === build groups exactly as Dict{String, Vector} ===
# (Vector element type may vary; your plot function handles Vector or SimulationState)
groups = Dict{String, Vector}()

for (hint, label) in patterns
    matches = glob("*$(hint)*final*.bin", data_dir)
    isempty(matches) && (println("No matches for $label"); continue)

    selected = first(matches, min(length(matches), max_files))

    loaded = Any[]  # Vector (concrete: Vector{Any}) — compatible with your method sig
    for f in selected
        try
            st = load_state(f)  # SimulationState
            push!(loaded, st)
        catch e
            @warn "Failed to load $(basename(f))" error=e
        end
    end

    groups[label] = loaded
    println("Loaded $(length(loaded)) runs for $label")
end


plt_cdf = plot_cdf_longest_by_type(groups; λ=0.6, mode=:geodesic, use_ccdf=false)
savefig(plt_cdf, joinpath(out_dir, "longest_strand_cdf.png"))

println("✅ Plots written to $out_dir")
