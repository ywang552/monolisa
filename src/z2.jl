using Random, Plots, Dates


# ----- output + config (your style) -----
out_dir = joinpath(pwd(), "plots", "tmp"); mkpath(out_dir)
config = Config(
    ft=0.0005, box_size=10, nn_restriction=3, box_capacity=10,
    monomer_radius=1.0, grid_size=2000, max_monomers=50_000,
    file_path="Claudins/wt2_newsep.txt", grid_overlay=false
)
name_noext, _ = splitext(basename(config.file_path))
stamp = Dates.format(now(), "yyyymmdd_HHMMSS")
seeds = [101, 202, 303, 404, 505]

# ----- CDF overlay plot using your helper -----
plt_cdf = plot(title="Face-size CDFs", xlabel="Face area", ylabel="CDF", legend=:bottomright)

for s in seeds
    Random.seed!(s)
    state = run(config)
    backbones = compute_backbone(state; λ=0.6, mode=:geodesic)

    # strand plot per trial (your prefix/args)
    prefix = joinpath(out_dir, name_noext * "_seed$(s)_")
    generate_plots(state, config; bbs = backbones, output_prefix = prefix*"arcs_$(stamp)", show_contour=true, tm_style=:nothing)

    # face areas + CDF overlay using your functions
    res = Faces.compute_enclosed_faces(state; area_floor=5.0, drop_outer=true)
    plot_face_area_cdf(res; logx=true, p=plt_cdf, min_area=5.0, label="seed=$(s)")
end

overlay_path = joinpath(out_dir, name_noext * "_overlay_faceCDF_$(stamp).png")
savefig(plt_cdf, overlay_path)
@info "Done — CDF saved at $overlay_path; strand images saved with prefix $(joinpath(out_dir, name_noext))*"
