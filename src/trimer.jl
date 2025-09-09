using Plots
include("plot_results.jl")   # load your plotting utilities

# --- define a simple trimer ---
x = [0.0, .4, .20]              # x-coordinates
y = [0.0, 0, .4*sqrt(3)/2]              # y-coordinates
rot = [0, 273, 18]               # spin angles (degrees)
r = 0.2                          # monomer radius

# edges connecting monomers (M1-M2-M3 chain)
edges = [(1,2), (2,3)]


# --- plot ---
plt = plot_monomers_lod(
    x, y;
    rotation=rot,
    monomer_radius=r,
    show_grid=false,
    show_orientation=true,
    lod=:detail,
    draw_contacts=true
)

# overlay backbone edges
draw_backbone!(plt, x, y, edges; lc=:black, lw=2)

display(plt)
savefig("trimer_double_contact.png")