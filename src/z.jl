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
    "Claudins/c4_7kp4.txt",
    "Claudins/hc15AF.txt",
    "Claudins\\wt2_newsep.txt",
    "Claudins\\hc5AF.txt"
]

# --- seeds ---
seeds = 101:101:505



# --- main loop ---
function main()
        # pick a palette you like (extend if you have more files)
    palette = [:red, :blue, :green, :purple, :orange, :teal, :brown, :magenta]

    # consistent color per file_path
    color_map = Dict{String, Symbol}()
    for (i, fp) in enumerate(claudin_files)
        color_map[fp] = palette[mod1(i, length(palette))]
    end

    # make the base CDF plot once
    plt_cdf = plot(
        xlabel = "Face area",
        ylabel = "CDF",
        title = "Face-size CDFs",
        legend = :bottomright,
    )
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
                max_monomers=50_000,
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
            col = color_map[file_path]

            plot_face_area_cdf(res; c = col, logx=true, p=plt_cdf, min_area=5.0, label="$(cn) seed=$(seed)")
            display(plt_cdf)
        end
    end
end 

main()