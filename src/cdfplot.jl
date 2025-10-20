using Serialization

folder = "data"
keywords = ["wt2","hc15"]

# 1. Get all serialized files (you can narrow this down by extension if needed)
all_files = readdir(folder; join=true)


# color_map = Dict{String, Symbol}()
# for (i, fp) in enumerate(keywords)
#     color_map[fp] = palette[mod1(i, length(palette))]
# end



plt_cdf = plot(
    xlabel = "Mesh area (nm²)",
    ylabel = "CDF",
    title = "Mesh-size CDFs",
)

palette = [:red, :blue, :green, :purple, :orange, :teal, :brown, :magenta]



for k in keywords 
    matched_files = filter(f -> occursin(k, lowercase(basename(f))), all_files)
    println("Matched files:")
    for f in matched_files
        println(f)
    end
    loaded_states = [deserialize(f) for f in matched_files]
    println("✅ Loaded $(length(loaded_states)) serialized states.")
    clor = color_map[k]
    for s in eachindex(loaded_states)[1]
        # lb = "wt2 #$(s)"
        lb = k*"$(s)"

        backbones = compute_backbone(loaded_states[s]; λ=0.6, mode=:geodesic)
        generate_plots(loaded_states[s], config; bbs = backbones, output_prefix = k*"$(s)", show_contour=true, tm_style=:nothing)

        # res = Faces.compute_enclosed_faces(loaded_states[s]; area_floor=5., drop_outer=true)
        # plot_face_area_cdf(res; c = clor, logx=true, p=plt_cdf, min_area=5., label=lb)
    end 
end 
plot!(plt_cdf; ylim=(0.7, 1.0))
plt_cdf


plot!(plt_cdf; ylim=(0.0, 0.7))
plt_cdf

plot!(plt_cdf; ylim=(0.0, 1.0))
plt_cdf

# savefig("asd.png")


