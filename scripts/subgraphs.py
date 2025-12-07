import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import re
from pathlib import Path

# read PPI network
ppi_path = "data/consensus_ppi_bioplex_biogrid_intact_huri_edgelist.tsv"

edges = pd.read_csv(
    ppi_path,
    sep="\t",
    header=None,
    usecols=[0, 1],
    names=["u", "v"],
    dtype=str
)
# remove self-loops
edges = edges[edges["u"] != edges["v"]]

# read Localization Data
loc_df = pd.read_csv("data/location_with_uniprot.tsv", sep="\t", dtype=str)

# get all unique proteins
#proteins = loc_df["uniprot_id"].unique().tolist()
#locations = loc_df["location"].unique().tolist()

#print(proteins)
#print(locations)

#Mapping dictionaries
loc_to_proteins = (         # Proteins in a single location
    loc_df
    .groupby("location")["uniprot_id"]
    .apply(set)
    .to_dict()
)



protein_to_locations = (    # Locations of a single protein
    loc_df
    .groupby("uniprot_id")["location"]
    .apply(lambda s: sorted(set(s)))
    .to_dict()
)

#print(loc_to_proteins["Intermediate filaments"])
#print(protein_to_locations["Q9Y2L8"])

#NetworkX Graph

G = nx.Graph()
G = nx.from_pandas_edgelist(edges, "u", "v")

all_nodes = set(G.nodes()) #all unique nodes in the PPI network

# PPI nodes without localisation
nodes_without_loc = all_nodes.difference(loc_df["uniprot_id"])
#print(nodes_without_loc)
Path("ppi_nodes_without_loc.tsv").write_text("\n".join(sorted(nodes_without_loc)))


location_subgraphs = {}

for loc, prots in loc_to_proteins.items(): #Iterates through all locations and their proteins
    #only proteins which are in the PPI network AND localization data
    nodes_here = list(all_nodes.intersection(prots))

    if not nodes_here:
        print(f"no nodes for {loc} location in the PPI")
        continue

    H_loc = G.subgraph(nodes_here).copy()   #extract subgraph via intersecting nodes
    location_subgraphs[loc] = H_loc         #save subgraph in dictionary with location as key and its subgraph as item

    print(f"{loc} Subnetwork: {H_loc.number_of_nodes()} nodes, {H_loc.number_of_edges()} edges")

#returns H_loc dictionary (Location: Subgraph)





#Check if node is in Subnetwork
print("P05981" in location_subgraphs["Plasma membrane"].nodes())
print("P05981" in location_subgraphs["Endoplasmic reticulum"].nodes())  
#should both be true as P05981 is multilocalized

#Search Edges of a  specific Protein in the network
protein = "Q9UJQ7"
edges_with_protein = list(G.edges(protein))

print(f"Protein {protein} Edges in Network:")
for e in edges_with_protein:
    print(e)

'''
#2D Plot Subgraphs
MAX_NODES = 50 

for loc, H in location_subgraphs.items():
    if H.number_of_nodes() > MAX_NODES:
        deg = dict(H.degree())
        top_nodes = [
            n for n, _ in sorted(deg.items(), key=lambda kv: kv[1], reverse=True)[:MAX_NODES]
        ]
        H_plot = H.subgraph(top_nodes).copy()
    else:
        H_plot = H

    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(H_plot, seed=42)

    degH = dict(H_plot.degree())
    sizes = np.array([max(d, 1) for d in degH.values()], dtype=float)
    sizes = 80 * np.sqrt(sizes)

    with_labels = H_plot.number_of_nodes() <= 60

    nx.draw_networkx(
        H_plot,
        pos=pos,
        with_labels=with_labels,
        node_size=sizes,
        font_size=8 if with_labels else 0,
        width=0.6,
    )

    plt.axis("off")
    plt.title(
        f"{loc} subgraph: {H_plot.number_of_nodes()} nodes, "
        f"{H_plot.number_of_edges()} edges"
    )

    out_png = Path(f"/Plots/ppi_{loc.replace(' ', '_')}.png")
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved plot for {loc} to {out_png}")

'''

from geometry import sample_ellipsoid #gets sample ellipsoid function from geometry.py

print(sample_ellipsoid((0.0, 0.0, 0.0), (4.0, 2.0, 2.0), n_points=1, rng=42)[0])

#create pos3d
def layout_organelle_random(H_loc, center, axes, rng=42):
    nodes = list(H_loc.nodes())
    points = sample_ellipsoid(center, axes, n_points=len(nodes), rng=rng) #sample points from geometry function
    # mapping: protein to 3D coordinates
    pos3d = {n: p for n, p in zip(nodes, points)} #build pos3d dictionary -> {protein_id: np.array([x, y, z]), ...}
    return pos3d

#Example with Nuclear membrane
center_nm = (0.0, 0.0, 0.0)
axes_nm = (4.0, 2.0, 2.0)

H_nm = location_subgraphs["Nuclear membrane"]
pos3d_nm = layout_organelle_random(H_nm, center_nm, axes_nm)

print(pos3d_nm)

#Plot in 3D

import plotly.graph_objects as go
import numpy as np

def plot_3d_network(H, pos3d, title="3D network"):
    # 3D Layout with NetworkX
    
    #pos3d is an dictionary with all nodes as keys and an (x,y,z) coordinate array as item
    #pos3d = {protein_id: np.array([x, y, z]), ...} -> pos3d[node] = (x, y, z)
    
    #pos3d = nx.spring_layout(H, dim=3, seed=42) -> this function should be replaced with our placing logic (geometry.py)
                                                
    #print(pos3d)

    # Edges
    edge_x = []
    edge_y = []
    edge_z = []

    for u, v in H.edges():
        x0, y0, z0 = pos3d[u]
        x1, y1, z1 = pos3d[v]
        edge_x += [x0, x1, None]
        edge_y += [y0, y1, None]
        edge_z += [z0, z1, None]

    edge_trace = go.Scatter3d(
        x=edge_x, y=edge_y, z=edge_z,
        mode="lines",
        line=dict(width=1),
        hoverinfo="none",  #Hover info of Edges
    )

    # Nodes
    node_x = []
    node_y = []
    node_z = []
    texts = []
    degrees = []

    for n in H.nodes():
        x, y, z = pos3d[n]
        node_x.append(x)
        node_y.append(y)
        node_z.append(z)
        texts.append(n)             # Hover-Text = Protein-ID
        degrees.append(H.degree[n])

    node_sizes = np.array(degrees, dtype=float)
    node_sizes = 4 + 5 * np.sqrt(np.maximum(node_sizes, 1))

    node_trace = go.Scatter3d(
        x=node_x, y=node_y, z=node_z,
        mode="markers",
        marker=dict(
            size=node_sizes,
            opacity=0.8,
        ),
        hovertext=texts,
        hoverinfo="text",
    )

    # build shape
    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(
        title=title,
        showlegend=False,
        margin=dict(l=0, r=0, b=0, t=40),
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
        )
    )

    fig.show()

H = location_subgraphs["Nuclear membrane"]

'''
MAX_NODES = 50

if H.number_of_nodes() > MAX_NODES:
    deg = dict(H.degree())
    top_nodes = [
        n for n, _ in sorted(deg.items(), key=lambda kv: kv[1], reverse=True)[:MAX_NODES]
    ]
    H = H.subgraph(top_nodes).copy()
'''

plot_3d_network(
    H,
    pos3d_nm,
    title=f"Endoplasmic reticulum (3D): {H.number_of_nodes()} nodes, {H.number_of_edges()} edges"
)