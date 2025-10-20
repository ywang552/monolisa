using Serialization, Plots
using Plots
using StatsBase          # for Histogram
# For box/violin plots:
# using StatsPlots         # pkg: StatsPlots.jl

const TRUE_COLORS = Dict("C2"=>:green, "C4"=>:red, "C5"=>:blue, "C15"=>:orange)

# --- infer claudin type from filename (customize patterns as needed) ---
function infer_type(fname::AbstractString; patterns=Dict(
        "hc5AF" => "C5",
        "wt2_newsep" => "C2",
        "hc15" => "C15",
        "c4_7" => "C4",
    ))
    b = basename(fname)
    # 1) direct match on substrings (robust & simple)
    for (needle, label) in patterns
        if occursin(needle, b)
            return label
        end
    end
    # 2) fallback: try regex like _C<number>_ or .C<number>.
    if (m = match(r"[_\.](C\d+)[_\.]", b)) !== nothing
        return m.captures[1]
    end
    return "Other"
end

# --- load all .bin states from a folder and group by inferred type ---
function load_states_grouped(dir::AbstractString; patterns=Dict{String,String}())
    groups = Dict{String,Vector}()
    for f in readdir(dir; join=true)
        endswith(f, ".bin") || continue
        cl = infer_type(f; patterns=patterns)
        st = deserialize(f)
        push!(get!(groups, cl, Vector()), st)
    end
    return groups
end

# --- your area extractor (uses Faces) ---
get_areas(state; area_floor=5.0) = begin
    res = Faces.compute_enclosed_faces(state;
        area_floor=area_floor, drop_outer=true,
        normalize_orientation=true, return_abs=true)
    res.areas
end

# --- ECDF on shared log grid + cluster plot (from previous message) ---
using StatsBase
function cdf_on_log_grid(areas::Vector{Float64}, xgrid_log::AbstractVector{Float64})
    logA = log10.(areas)
    ec   = ecdf(logA)
    @inbounds [ec(x) for x in xgrid_log]
end

function plot_cdf_clusters_by_type(data_by_type::Dict{String,Vector};
        colors=Dict{String,Symbol}(), area_floor=8.0,
        use_ccdf::Bool=true, focus_top::Union{Nothing,Tuple{Float64,Float64}}=(0.7,1.0))

    # collect all areas for shared grid
    all_logA = Float64[]
    pertype_areas = Dict{String,Vector{Vector{Float64}}}()
    for (cl, runs) in data_by_type
        per = Vector{Vector{Float64}}()
        for st in runs
            A = isa(st, Vector) ? st : get_areas(st; area_floor=area_floor)
            A = filter(>=(area_floor), A)
            isempty(A) && continue
            push!(per, A)
            append!(all_logA, log10.(A))
        end
        pertype_areas[cl] = per
    end
    isempty(all_logA) && error("No faces after filtering area_floor=$(area_floor).")

    # shared x-grid (log-space)
    xmin, xmax = extrema(all_logA)
    xgrid_log  = range(xmin, xmax; length=500)
    xgrid_area = 10 .^ xgrid_log .* 1e-6

    plt = plot(; xlabel="Mesh area (μm²)", ylabel=use_ccdf ? "CDF" : "CDF",
               xscale=:log10, legend=:bottomright, grid=true,
               title="Mesh-size CDFs by Claudin")

    for (cl, runs) in sort(collect(pertype_areas); by=first)
        col = get(colors, cl, get(Dict("C4"=>:red,"C2"=>:green,"C5"=>:blue,"C15"=>:orange), cl, :purple))

        Ys = Float64[]
        counts = Int[]
        for A in runs
            y = cdf_on_log_grid(A, xgrid_log)
            if use_ccdf; y = y; end
            plot!(plt, xgrid_area, y; color=col, alpha=0.25, lw=1, label = false)
            append!(Ys, y); push!(counts, length(y))
        end
        isempty(counts) && continue

        m = length(xgrid_log)
        Y = reshape(Ys, m, :)
        μ = mapslices(mean, Y; dims=2)[:]
        σ = mapslices(std,  Y; dims=2)[:]
        plot!(plt, xgrid_area, μ; color=col, lw=3, label=cl, ribbon=σ, fillalpha=0.15)
    end
            plot!(plt, xgrid_area, y; color=col, alpha=0.25, lw=1, label = false)

    if focus_top !== nothing
        ylims!(plt, focus_top...)
    else
        ylims!(plt, 0, 1)
    end
    return plt
end

# Collect all areas (nm²) for one claudin label by scanning files and computing faces
function collect_areas_for_claudin(data_dir::AbstractString, patterns::Dict{String,String},
                                   target::String; area_floor::Float64=8.0)
    areas = Float64[]
    for f in readdir(data_dir; join=true)
        endswith(f, ".bin") || continue
        cl = infer_type(f; patterns=patterns)   # you already have this
        cl == target || continue
        st = deserialize(f)
        # compute areas from faces (keep it simple; drop outer; make areas positive)
        res = Faces.compute_enclosed_faces(st; area_floor=area_floor, drop_outer=true, return_abs=true)
        append!(areas, res.areas)
    end
    filter!(x -> x >= area_floor, areas)
    areas .= areas.*1e-6
    return areas
end
# Make one histogram per claudin; returns a Dict of plots
function plot_histograms_by_claudin(data_dir::AbstractString, patterns::Dict{String,String};
                                    area_floor::Float64=8.0, nbins::Int=60, density::Bool=true)
    labels = ["C2","C4","C5","C15"]
    plots  = Dict{String,Plots.Plot}()

    for lab in labels
        A = collect_areas_for_claudin(data_dir, patterns, lab; area_floor=area_floor)
        s = round(area_floor*1e-6, digits=6)


        p = histogram(A.*1e-6;
            bins=nbins,
            normalize = density,          # <- Bool, not :pdf / :none
            xlabel="Mesh area (μm²)",
            ylabel = density ? "Density" : "Count",
            title=lab*"_N = $(length(A))_Δ > $(s) μm²",
            legend=false,
            grid=true,
            yscale = :log10)
        plots[lab] = p
    end
    return plots
end

# Optional: assemble a 2×2 grid of the four histograms
function plot_histograms_grid(data_dir::AbstractString, patterns::Dict{String,String};
                              area_floor::Float64=8.0, nbins::Int=60, density::Bool=true)
    ps = plot_histograms_by_claudin(data_dir, patterns; area_floor=area_floor, nbins=nbins, density=density)
    return plot(ps["C2"], ps["C4"], ps["C5"], ps["C15"]; layout=(2,2), size=(1000,700))
end


function collect_areas_for_claudin(data_dir::AbstractString, patterns::Dict{String,String},
                                   target::String; area_floor::Float64=8.0)
    areas = Float64[]
    for f in readdir(data_dir; join=true)
        endswith(f, ".bin") || continue
        cl = infer_type(f; patterns=patterns)
        cl == target || continue
        st = deserialize(f)
        res = Faces.compute_enclosed_faces(st; area_floor=area_floor, drop_outer=true, return_abs=true)
        append!(areas, res.areas)
    end
    filter!(x -> x >= area_floor, areas)
    return areas
end

function plot_box_whisker_by_claudin(data_dir::AbstractString, patterns::Dict{String,String};
                                     area_floor::Float64=8.0)
    labels = ["C2","C4","C5","C15"]
    xs, ys, cols = String[], Vector{Vector{Float64}}(), Symbol[]

    for lab in labels
        A = collect_areas_for_claudin(data_dir, patterns, lab; area_floor=area_floor)
        isempty(A) && continue
        push!(xs, lab)
        push!(ys, A)
        push!(cols, get(TRUE_COLORS, lab, :gray))
    end

    # One series per claudin; set all visuals to :match so each series is a single color
    plt = boxplot(xs, ys;
        color=cols, fillalpha=0.35, legend=false, linewidth=1.5,
        whiskercolor=:match,     # whiskers same as box color
        outliercolor=:match,     # outlier dots same as box color
        mediancolor=:black,      # median line stands out
        xlabel="Claudin", ylabel="Area", title="Mesh-area spread by claudin", grid=true
    )
    return plt
end



# ------------------ RUN IT ------------------
# Folder with your .bin states
data_dir = "data"

# Customize how filenames map to claudin labels, if needed
patterns = Dict(
    # "C4" => "C4",
    # "C2" => "C2",
    "hc5AF" => "C5",
    "wt2_newsep" => "C2",
    "hc15" => "C15",
    "c4_7" => "C4",

    # add more keys if filenames use other hints, e.g. "claudin4" => "C4"
)
area_floor  = 8.
groups = load_states_grouped(data_dir; patterns=patterns)
# pgrid = plot_histograms_grid(data_dir, patterns; area_floor=10., nbins=60, density=false)
ps = plot_histograms_by_claudin(data_dir, patterns, density = false, area_floor = area_floor)


for p in keys(ps)
    savefig(ps[p], "figs\\cdf\\$(p)_frequency_$(area_floor).png")
end 

# patterns = Dict("hc5AF"=>"C5","wt2_newsep"=>"C2","hc15"=>"C15","c4_7"=>"C4")


labels = ["C2","C4","C5","C15"]
xs, ys, cols = String[], Vector{Vector{Float64}}(), Symbol[]

for lab in labels
    A = collect_areas_for_claudin(data_dir, patterns, lab; area_floor=100.)
    isempty(A) && continue
    push!(xs, lab)
    push!(ys, A)
    push!(cols, get(TRUE_COLORS, lab, :gray))
end

# using DataFrames, CSV
# p = [fill("C2", length(ys[1])), fill("C4", length(ys[2])), fill("C5", length(ys[3])), fill("C15", length(ys[4]))]
# k = vcat(ys...)
# k = k.*1e-6
# # DataFrame(Claudin = )
# df = DataFrame(Claudin = vcat(p...),
#                Area    = k)
# CSV.write("figs/whiskerbox/areas.csv", df)



# using PlotlyJS
# data = GenericTrace[]
# for i in 1:length(labels)
#     trace = box(;y=ys[i],
#                     name=xs[i],
#                     boxpoints="all",
#                     jitter=0.5,
#                     whiskerwidth=0.2,
#                     fillcolor="cls",
#                     marker_size=2,
#                     line_width=1)
#     push!(data, trace)
# end

# t = "s"
# layout = Layout(;title=t,
#                     yaxis=attr(autorange=true, showgrid=true, zeroline=true,
#                             dtick=5, gridcolor="rgb(255, 255, 255)",
#                             gridwidth=1,
#                             zerolinecolor="rgb(255, 255, 255)",
#                             zerolinewidth=2),
#                     margin=attr(l=40, r=30, b=80, t=100),
#                     paper_bgcolor="rgb(243, 243, 243)",
#                     plot_bgcolor="rgb(243, 243, 243)",
#                     showlegend=false)
# PlotlyJS.plot(data, layout)

# for v in values(patterns)
#     savefig(ps[v], "figs\\cdf\\$(v)_histogram_density.png")
# end 

# plt = plot_cdf_clusters_by_type(groups;
#     colors=Dict("C4"=>:red, "C2"=>:green),
#     area_floor=8.2,
#     use_ccdf=true,
#     focus_top=(0.7, 1.0),
# )[1]

# display(plt)
# # savefig(plt, "figs/cdf/Focus Top.png")


# plt = plot_cdf_clusters_by_type(groups;
#     colors=Dict("C4"=>:red, "C2"=>:green),
#     area_floor=8.2,
#     use_ccdf=true,
#     focus_top=(0.0, 1.0), 
# )[1]

# plot_hist_by_type(groups)
# savefig(plt, "figs/cdf/Full.png")
