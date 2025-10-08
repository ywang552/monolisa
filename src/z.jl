using Serialization, Dates

function make_state_filename(state, config; seed::Int, prefix="data")
    name = splitext(basename(config.file_path))[1]
    N = state.MN
    tag = "mono"
    ts = Dates.format(now(), "yyyymmdd_HHMM")
    fname = "$(tag)_$(name)_N$(N)_seed$(seed)_$(ts).bin"
    return joinpath(prefix, fname)
end


# --- input W matrices ---
claudin_files = [
    "Claudins\\wt2_newsep.txt",
    "Claudins\\hc5AF.txt",
    "Claudins/hc15AF.txt",
    "Claudins/c4_7kp4.txt",
]

# --- seeds ---
seeds = 1:1



plt_cdf = plot(title="Face-size CDFs", xlabel="Face area", ylabel="CDF", legend=:bottomright)

# --- main loop ---
function main()
    for file_path in claudin_files
        for seed in seeds
            println("\n=== Running $(basename(file_path)) | seed=$(seed) ===")
            println(file_path)
            config123 = Config(
                ft=0.0005,
                box_size=10,
                nn_restriction=3,
                box_capacity=10,
                monomer_radius=1,
                grid_size=2000,
                max_monomers=3000,
                # file_path=ARGS[1],
                file_path=file_path,
                grid_overlay = false,
                prog = 1
            )
            Random.seed!(seed)
            state = F(config123)

            save_name = make_state_filename(state, config123, seed = seed)
            serialize(save_name, state)
            println("✅ Saved: $save_name  |  MN=$(state.MN)")
            res = Faces.compute_enclosed_faces(state; area_floor=5.0, drop_outer=true)

            cn = replace(splitdir(file_path)[end], ".txt" => "")

            plot_face_area_cdf(res; logx=true, p=plt_cdf, min_area=5.0, label="seed=$(cn) $(seed)")
            display(plt_cdf)
        end
    end
end 

main()