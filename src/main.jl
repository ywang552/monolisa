# ==============================================
# File: main.jl
# Version: v1.0.0
# Description: Entry point for MonoLisa simulation.
# Author: [Youchuan Wang]
# Date: [11/17/2024]
# ==============================================

include("dependencies.jl")
include("initialization.jl")
include("simulation.jl")
include("utils.jl")
include("plot_results.jl")
include("structure_analysis.jl")



"""
Sets up the necessary folders for the project.
Creates folders if they do not already exist.
"""
function setup_project_folders()
    folders = ["saved_states", "saved_states\\minimal_states", "saved_states\\full_states", "Claudins", "plots", "logs"]
    for folder in folders
        folder = pwd()*"\\"*folder
        if !isdir(folder)
            mkdir(folder)
            println("Created folder: $folder")
        else
            println("Folder already exists: $folder")
        end
    end
end

# setup_project_folders()

config = Config(
    ft=0.0005,
    box_size=10,
    nn_restriction=3,
    box_capacity=10,
    monomer_radius=1,
    grid_size=2000,
    max_monomers=10000,
    # max_monomers=100,
    # file_path=ARGS[1],
    file_path="Claudins/8Cpop_grids.txt",
    grid_overlay = false,
    prog = 1
)

function main()
    println("starting the placement...")
    # Pass configuration to simulation
    state = run(config)

    # Save results or generate plots
    return state
end


"""
    report_overconnected(state; max_deg=3, use_geometric=false, contact_scale=1.02)

Check each monomer's number of neighbors and report any with degree > max_deg.
"""
function report_overconnected(state;
    max_deg::Int=2,
    use_geometric::Bool=false,
    contact_scale::Real=1.02
)

    x, y, r = state.x_coords, state.y_coords, state.radius
    xmin, xmax = minimum(x), maximum(x)
    ymin, ymax = minimum(y), maximum(y)
    pad = 2r
    xlim = (xmin - pad, xmax + pad)
    ylim = (ymin - pad, ymax + pad)
    edge_list = find_contact_edges(x, y, r; xlim=xlim, ylim=ylim,
                                    contact_scale=contact_scale)


    # --- degree map ---
    deg = Dict{Int,Int}()
    for (u,v) in edge_list
        u == v && continue
        deg[u] = get(deg,u,0) + 1
        deg[v] = get(deg,v,0) + 1
    end

    offenders = [i for (i,d) in deg if d > max_deg]

    println("Neighbor-degree check (threshold = $max_deg)")
    println("Total monomers: ", length(state.x_coords),
            " | edges considered: ", length(edge_list))
    if isempty(offenders)
        println("✓ No monomers exceed degree $max_deg.")
    else
        println("⚠ ", length(offenders), " monomer(s) exceed degree $max_deg:")
        for i in offenders
            println("  monomer $i → degree $(deg[i])")
        end
    end

    return (; max_deg, offenders, degrees=deg, edge_list)
end




safe_stamp() = Dates.format(now(), "yyyymmdd-HHMMSS-sss")
stamp = safe_stamp()

state = main()


ms = create_minimal_state(state)
ms.edges
state.edges
ms.junctions
ms.endpoints

# print(length(state.x_coords))
out_dir = joinpath(pwd(), "plots", "tmp")
mkpath(out_dir)

name_noext, _ = splitext(basename(config.file_path))
prefix = joinpath(out_dir, name_noext*"_123")
overlay_path = prefix * "overlay_curvature.png"

generate_plots(state, config; output_prefix = prefix*"arcs_$(stamp)", show_contour=true, tm_style=:nothing)

# import DelaunayTriangulation as DT

# # Shoelace area for a polygon given as Vector{Tuple{<:Real,<:Real}}.
# shoelace(coords) = begin
#     n = length(coords); n < 3 && return 0.0
#     s = 0.0
#     @inbounds for i in 1:n
#         xi, yi = coords[i]
#         xj, yj = coords[i == n ? 1 : i+1]
#         s += xi*yj - xj*yi
#     end
#     0.5*abs(s)
# end

# function voronoi_total_area(state; padding=0.05)
#     xs, ys = state.x_coords, state.y_coords
#     tri = DT.triangulate([(x,y) for (x,y) in zip(xs, ys)])

#     xmin, xmax = extrema(xs); ymin, ymax = extrema(ys)
#     w = max(xmax - xmin, eps()); h = max(ymax - ymin, eps())
#     δx, δy = w*padding, h*padding
#     xlims = (xmin - δx, xmax + δx)
#     ylims = (ymin - δy, ymax + δy)

#     vorn = DT.voronoi(tri; xlims, ylims, clip=true)
#     bbox = (xlims[1], xlims[2], ylims[1], ylims[2])

#     total = 0.0
#     for (poly_idx, _) in DT.get_polygons(vorn)
#         coords = DT.get_polygon_coordinates(vorn, poly_idx, bbox)
#         total += shoelace(coords)
#     end
#     bbox_area = (xlims[2]-xlims[1]) * (ylims[2]-ylims[1])
#     total, bbox_area, total / bbox_area
# end

# sum_area, clip_area, ratio = voronoi_total_area(state)
# @info "Voronoi sum = $sum_area, clip box = $clip_area, ratio = $ratio"


# # Andrew’s monotone chain (no extra deps)
# function convex_hull(points::Vector{Tuple{T,T}}) where {T<:Real}
#     pts = sort(points)
#     cross(o,a,b) = (a[1]-o[1])*(b[2]-o[2]) - (a[2]-o[2])*(b[1]-o[1])
#     lower = Tuple{T,T}[]
#     for p in pts
#         while length(lower) >= 2 && cross(lower[end-1], lower[end], p) <= 0
#             pop!(lower)
#         end
#         push!(lower, p)
#     end
#     upper = Tuple{T,T}[]
#     for p in Iterators.reverse(pts)
#         while length(upper) >= 2 && cross(upper[end-1], upper[end], p) <= 0
#             pop!(upper)
#         end
#         push!(upper, p)
#     end
#     vcat(lower[1:end-1], upper[1:end-1])
# end

# pts = [(x,y) for (x,y) in zip(state.x_coords, state.y_coords)]
# hull = convex_hull(pts)
# hull_area = shoelace(hull)
# @info "Convex hull area = $hull_area"


# plt = plot_voronoi_edges_DT(state)
# display(plt)
# savefig(plt,"xd.png")
# # Fast path: uses state.edges from the simulation
# res = report_overconnected(state; max_deg=3)

# # Cross-check against geometric contacts (same logic as backbone plotting)
# res_geom = report_overconnected(state; max_deg=3, use_geometric=true, contact_scale=1.02)


x = state.x_coords
y = state.y_coords;
state.edges
# using CSV, DataFrames

# df_a = DataFrame(x = state.x_coords, y = state.y_coords, r = state.rotation)
# CSV.write("3D_info.csv", df_a)
