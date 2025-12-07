import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from pathlib import Path

import os
print("CWD =", os.getcwd())

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

whole_cell_G = nx.Graph()
whole_cell_G = nx.from_pandas_edgelist(edges, "u", "v")

all_nodes = set(whole_cell_G.nodes()) #all unique nodes in the PPI network

# PPI nodes without localisation
nodes_without_loc = all_nodes.difference(loc_df["uniprot_id"])
nodes_with_loc = all_nodes.intersection(loc_df["uniprot_id"])
#print(nodes_without_loc)
#Path("ppi_nodes_without_loc.tsv").write_text("\n".join(sorted(nodes_without_loc)))

#from geometry import sample_ellipsoid

#NetworkX Attributes e.g Graph.nodes[node]["attribute"] = "red"

for node in nodes_with_loc:
    whole_cell_G.nodes[node]["locations"] = {} #should be a dictionary

    for location in protein_to_locations[node]:
         #3D sampling needs to be replaced with our placing logic -> 3D space
         #fill dictionary with locations and 3D coordinates
         #G.nodes[node]["locations"][location] = sample_ellipsoid((0.0, 0.0, 0.0), (4.0, 2.0, 2.0), n_points=1, rng=42)[0]
         continue






if __name__ == "__main__":
    '''
    #Print-Loop for all nodes and its attributes of the Network

    for node, attrs in G.nodes(data=True):
        print(node, attrs)
    '''

    '''
    #Print-Loop for all nodes and its attributes of a Subnetwork "Cytosol"

    for node, attrs in location_subgraphs["Cytosol"].nodes(data=True):
        print(node, attrs)
    '''

    #edge list with locations

    edge_nodes_with_attrs = [
        ((u, dict(whole_cell_G.nodes[u])), (v, dict(whole_cell_G.nodes[v])))
        for u, v in whole_cell_G.edges()
    ]

    for pair in edge_nodes_with_attrs[:50]:
        #print(pair)
        continue

    #print(len(edge_nodes_with_attrs))

    #print(G.nodes["Q9BUF5"]["locations"]["Microtubules"]) #Acess coordinates of a Node/Protein
    #print(nx.get_node_attributes(G, "locations").keys())