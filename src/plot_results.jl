using Plots
using Statistics

"""
Fast plotting of monomer positions with coordinate transformations.
"""
function fast_plot_monomers(
    x_coords, y_coords, boxSize, plot_grid, cl;
    marker_size=0.05, save_path=nothing, hide_title=true, hide_labels=true, hide_axes=true
)
    # Offset and scaling
    x_offset = minimum(x_coords)
    y_offset = minimum(y_coords)
    x_coords .= x_coords .- x_offset
    y_coords .= y_coords .- y_offset
    x_coords .= x_coords .* 0.37
    y_coords .= y_coords .* 0.37

    # Plot monomers
    scatter(
        x_coords, y_coords,
        color=:green, marker=:circle, ms=marker_size,
        fillcolor=:green, linecolor=:green, aspect_ratio=1, legend=false, grid=false, dpi=600
    )

    # Handle optional title and labels
    if !hide_title
        plot!(title="Monomer Placement: $cl ($(length(x_coords)) Monomers)")
    end
    if !hide_labels
        plot!(xlabel="x (nanometers)", ylabel="y (nanometers)")
    end
    if hide_axes
        plot!(xticks=nothing, yticks=nothing, framestyle=:none)
    end

    # Optionally plot the grid
    if plot_grid
        x_max = maximum(x_coords)
        y_max = maximum(y_coords)
        for x in 0:boxSize*0.37:x_max
            plot!([x, x], [0, y_max], lw=0.5)  # Vertical lines
        end
        for y in 0:boxSize*0.37:y_max
            plot!([0, x_max], [y, y], lw=0.5)  # Horizontal lines
        end
    end

    # Save the plot if a save path is provided
    if save_path !== nothing
        savefig(save_path)
        println("Plot saved to: $save_path")
    end
end


"""
Plots a heatmap for a given matrix.
"""
function plot_heatmap(matrix, title=nothing; save_path=nothing, bg_color=:transparent)
    heatmap(matrix, 
            title=title, 
            xlabel="X Axis", 
            ylabel="Y Axis", 
            color=cgrad([:blue, :red]), # You can choose other color maps or create a custom one
            clims=(minimum(matrix), maximum(matrix)))  # Set color limits based on data range
    
    if save_path !== nothing
        savefig(save_path)
        println("Heatmap saved to: $save_path")
    end
end

"""
Main function to generate and save plots for monomer placement and heatmaps.

Arguments:
- `state`: The SimulationState struct containing simulation data.
- `config`: The Config struct with simulation configuration parameters.
- `output_prefix`: Prefix for saving the plots. Defaults to "results/monolisa".
"""
function generate_plots(state::AbstractState, config, output_prefix=nothing)
    # Extract common properties
    x_coords = state.x_coords
    y_coords = state.y_coords
    W = state.W
    K = state.K

    # Handle box size and grid overlay
    box_size = state.box_size
    overlay = config.grid_overlay

    # File information
    file_name = basename(config.file_path)
    num_monomers = length(x_coords)
    # Generate plots
    fast_plot_monomers(
        x_coords, y_coords, box_size, overlay, file_name, 
        save_path="$(output_prefix)_$(num_monomers)_$(file_name)_monomers_placement.png"
    )
    plot_heatmap(
        W, "in_$(file_name)_$(num_monomers)", 
        save_path="$(output_prefix)_$(num_monomers)_$(file_name)_monomers_in.png"
    )
    plot_heatmap(
        K, "out_$(file_name)_$(num_monomers)", 
        save_path="$(output_prefix)_$(num_monomers)_$(file_name)_monomers_out.png"
    )

    println("All plots saved with prefix: $output_prefix")
end


# Entry point to use this script independently
if abspath(PROGRAM_FILE) == @__FILE__
    println("Running plotting script...")
    # Add code to load saved state/config and call `generate_plots` if running directly
    # Example:
    # include("your_simulation_saving_script.jl")
    # state, config = load_simulation_results("results/simulation_state.jld2")
    # generate_plots(state, config)
end
