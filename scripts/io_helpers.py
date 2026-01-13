"""I/O helpers: load PPI network + localization table."""

import pandas as pd
import networkx as nx
from pathlib import Path

def load_ppi_graph(ppi_path: str, sep: str = "\t") -> nx.Graph:
    """Load an undirected PPI graph from a two-column edge list file."""

    edges = pd.read_csv(ppi_path, sep="\t", names=["u", "v"], dtype=str)
    
    # remove self-loops
    edges = edges[edges["u"] != edges["v"]]
    
    G = nx.from_pandas_edgelist(edges, "u", "v")
    return G

def load_localizations(loc_path: str, sep: str = "\t", protein_col: str = "uniprot_id", location_col: str = "location",):
    """Load localization table and return mapping dicts.
        Returns
        -------
        loc_to_proteins:
            location -> set(protein IDs)
        protein_to_locations:
            protein ID -> sorted list of unique locations
    """

    loc_df = pd.read_csv(loc_path, sep=sep, dtype=str)

    #Mapping dictionaries
    
    loc_to_proteins = (         # Proteins in a single location
        loc_df
        .groupby(location_col)[protein_col]
        .apply(set)
        .to_dict()
    )

    protein_to_locations = (    # Locations of a single protein
    loc_df
    .groupby(protein_col)[location_col]
    .apply(lambda s: sorted(set(s)))
    .to_dict()
    )

    return loc_to_proteins, protein_to_locations, loc_df

def init_locations_attribute(G: nx.Graph) -> None:
    """Ensure every node has a locations dict attribute."""
    for node in G.nodes:
        G.nodes[node]["locations"] = {} #should be a dictionary

def get_nodes_without_location(G: nx.Graph, loc_df):
    all_nodes = set(G.nodes())
    return all_nodes.difference(loc_df["uniprot_id"])

def output_nodes_without_location(G: nx.Graph, loc_df, out_path: str):
    """creates .tsv file with all nodes, where no localiazation data is available"""
    all_nodes = set(G.nodes())
    nodes_without_loc = all_nodes.difference(loc_df["uniprot_id"])
    Path(out_path).write_text("\n".join(sorted(nodes_without_loc)))

def get_nodes_with_location(G: nx.Graph, loc_df):
    all_nodes = set(G.nodes())
    return all_nodes.intersection(loc_df["uniprot_id"])

def print_network(G):
    """Print-Loop for all nodes and its attributes of the Network"""
    for node, attrs in G.nodes(data=True):
        print(node, attrs)