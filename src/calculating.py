# Complete Code: Voronoi Tessellation with Colored Cells Based on Area

import numpy as np
import networkx as nx
from scipy.spatial import Voronoi, ConvexHull, KDTree, voronoi_plot_2d
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# Step 1: Generate synthetic monomer positions for a test case
np.random.seed(42)  # For reproducibility
num_monomers = 2000  # Adjust for larger or smaller dataset
monomer_positions = np.random.uniform(0, 600, size=(num_monomers, 2))  # Random (x, y) coordinates

# Step 2: Build a graph using KDTree for spatial connections
G = nx.Graph()
tree = KDTree(monomer_positions)

for i, pos in enumerate(monomer_positions):
    distances, indices = tree.query(pos, k=min(6, len(monomer_positions)))  # Increased connections
    
    for j in range(1, len(indices)):  # Skip the first index (itself)
        if distances[j] < 20:  # Increased threshold to ensure connections form clusters
            G.add_edge(tuple(pos), tuple(monomer_positions[indices[j]]))

# Step 3: Classify nodes in the graph
branching_nodes = {node for node in G.nodes if len(G[node]) >= 3}

# Step 4: Compute the Voronoi tessellation
vor = Voronoi(monomer_positions)

# Step 5: Compute Voronoi cell areas (ignoring infinite regions)
voronoi_areas = []
valid_regions = []
for region in vor.regions:
    if not -1 in region and len(region) > 0:  # Ignore infinite regions
        polygon = np.array([vor.vertices[i] for i in region])
        if len(polygon) >= 3:
            hull = ConvexHull(polygon)
            voronoi_areas.append(hull.volume)
            valid_regions.append(region)

# Step 6: Plot Voronoi cell size distribution
plt.figure(figsize=(8, 6))
plt.hist(voronoi_areas, bins=50, color="skyblue", edgecolor="black")
plt.xlabel("Voronoi Cell Area")
plt.ylabel("Frequency")
plt.title("Distribution of Voronoi Cell Sizes in Monomer Network")
plt.show()

# Step 7: Overlay Voronoi Diagram with Colored Cells Based on Area
fig, ax = plt.subplots(figsize=(12, 12))

# Convert positions to a dictionary for NetworkX
pos_dict = {tuple(map(float, node)): tuple(map(float, node)) for node in G.nodes}

# Draw the graph
nx.draw(G, pos=pos_dict, with_labels=False, node_size=2, edge_color="gray", ax=ax)

# Highlight branching nodes in red
nx.draw_networkx_nodes(G, pos=pos_dict, nodelist=list(map(tuple, branching_nodes)), 
                       node_color="red", node_size=6, label="Branching Nodes", ax=ax)

# Color Voronoi regions based on cell area
min_area, max_area = min(voronoi_areas), max(voronoi_areas)
norm = mcolors.Normalize(vmin=min_area, vmax=max_area)
colormap = cm.get_cmap("plasma")

for region, area in zip(valid_regions, voronoi_areas):
    polygon = np.array([vor.vertices[i] for i in region])
    ax.fill(polygon[:, 0], polygon[:, 1], color=colormap(norm(area)), alpha=0.5)

# Add colorbar for reference
sm = cm.ScalarMappable(norm=norm, cmap="plasma")
cbar = plt.colorbar(sm, ax=ax, label="Voronoi Cell Size")

plt.legend()
plt.title("Monomers: Voronoi Tessellation with Colored Cells Based on Area")
plt.show()

# Step 8: Display Summary Statistics for Voronoi Cell Sizes
{
    "Total Voronoi Cells": len(voronoi_areas),
    "Mean Voronoi Cell Area": np.mean(voronoi_areas),
    "Min Voronoi Cell Area": np.min(voronoi_areas),
    "Max Voronoi Cell Area": np.max(voronoi_areas),
    "Voronoi Areas (Sampled)": voronoi_areas[:10]
}

# Re-run the Voronoi plot separately to ensure it's displayed

fig, ax = plt.subplots(figsize=(12, 12))

# Draw the original Voronoi diagram
voronoi_plot_2d(vor, ax=ax, show_vertices=False, line_colors="blue", line_width=0.5, line_alpha=0.6)

# Highlight monomer positions
ax.scatter(monomer_positions[:, 0], monomer_positions[:, 1], color="red", s=10, label="Monomers")

plt.legend()
plt.title("Voronoi Tessellation of Monomer Network")
plt.show()



