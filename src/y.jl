# ---- 1) Longest-backbone length per state ----
# Robust to compute_backbone returning either a Vector of paths, or a tuple (paths, ...).
longest_backbone_len(state; λ=0.6, mode=:geodesic) = begin
    bb = compute_backbone(state; λ=λ, mode=mode)
    paths = bb isa Tuple ? bb[1] : bb
    isempty(paths) && return 0
    maximum(length(p) for p in values(paths)) - 1
end

# ---- 2) Collect lengths by type (same layout as your area collector) ----
# data_by_type :: Dict{String, Vector}  where each element is either a state or a precomputed Int
function collect_lens_by_type(data_by_type::Dict{String,Vector};
                              λ=0.6, mode=:geodesic)
    pertype = Dict{String,Vector{Int}}()
    for (cl, runs) in data_by_type
        vals = Int[]
        for st in runs
            if st isa Integer
                push!(vals, st)
            else
                push!(vals, longest_backbone_len(st; λ=λ, mode=mode))
            end
        end
        pertype[cl] = vals
    end
    return pertype
end

# ---- 3) CDF of longest-backbone length by claudin ----
# use_ccdf=true will show 1 - CDF.
function plot_cdf_longest_by_type(data_by_type::Dict{String,Vector};
        colors=Dict{String,Symbol}(),
        λ::Float64=0.6, mode::Symbol=:geodesic,
        use_ccdf::Bool=false, focus_top::Union{Nothing,Tuple{Float64,Float64}}=nothing)

    pertype = collect_lens_by_type(data_by_type; λ=λ, mode=mode)
    allvals = reduce(vcat, values(pertype))
    isempty(allvals) && error("No backbone lengths to plot.")

    # integer grid
    xmin, xmax = extrema(allvals)
    xgrid = collect(xmin:xmax)

    plt = plot(; xlabel="Longest strand length (edges)",
                ylabel=use_ccdf ? "CCDF" : "CDF",
                legend=:bottomright, grid=true,
                title="Longest-strand CDFs by Claudin")

    default_cols = Dict("C4"=>:red,"C2"=>:green,"C5"=>:blue,"C15"=>:orange)

    for (cl, vals) in sort(collect(pertype); by=first)
        isempty(vals) && continue
        col = get(colors, cl, get(default_cols, cl, :purple))

        # Empirical CDF on shared x-grid
        counts = zeros(Float64, length(xgrid))
        for v in vals
            idx = clamp(v - xmin + 1, 1, length(xgrid))
            counts[idx] += 1
        end
        cdf = cumsum(counts) ./ length(vals)
        y = use_ccdf ? (1 .- cdf) : cdf

        # per-run traces (faint) are optional; we keep it clean here.
        plot!(plt, xgrid, y; color=col, lw=3, label=cl)
    end

    if focus_top !== nothing
        ylims!(plt, focus_top...)
    else
        ylims!(plt, 0, 1)
    end

    return plt
end

# ---- 4) Histograms of longest-backbone length by claudin ----
# density=true normalizes to probability mass; bins default to one per integer.
function plot_hist_longest_by_type(data_by_type::Dict{String,Vector};
        colors=Dict{String,Symbol}(),
        λ::Float64=0.6, mode::Symbol=:geodesic,
        density::Bool=true)

    pertype = collect_lens_by_type(data_by_type; λ=λ, mode=mode)
    plots = Dict{String,Plots.Plot}()
    default_cols = Dict("C4"=>:red,"C2"=>:green,"C5"=>:blue,"C15"=>:orange)

    for cl in sort(collect(keys(pertype)))
        vals = pertype[cl]
        isempty(vals) && continue
        col = get(colors, cl, get(default_cols, cl, :purple))

        xmin, xmax = extrema(vals)
        # one bin per integer length
        nbins = max(1, xmax - xmin + 1)
        p = histogram(vals;
            bins=xmin:(xmax+1),         # integer edges → clean bin boundaries
            normalize=density,          # Bool
            xlabel="Longest strand length (edges)",
            ylabel = density ? "Density" : "Count",
            title = "$(cl)  N=$(length(vals))",
            legend=false, grid=true,
            linecolor=:auto, fillalpha=0.6, color=col)
        plots[cl] = p
    end
    return plots
end
