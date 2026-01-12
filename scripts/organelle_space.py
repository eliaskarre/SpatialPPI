import math
import random
from typing import List, Dict, Tuple, Optional, Callable
import numpy as np

import networkx as nx

import plotly.graph_objects as go

from geometry import *
from constraint_layout_classes import *

class Location:
    """
    Location (Organell) Object.

    Attributes:
    - name:             "Nucleus", "Cytosol" ..
    - nodes:            List of Proteins (Node-IDs) in this location
    - min_degree:       minimum number of degree per node
    - center:           3D coordinates of center (x, y, z)
    - constraint:       3D layout constrained by a certain shape
    """

    def __init__(
        self,
        name: str,
        node_ids: List[str],
        min_degree: int,
        center: Tuple[float, float, float],
        constraint: object,
        repulsion_strength = 1
    ):
        self.name = name
        self.node_ids = node_ids
        self.min_degree = min_degree
        self.center = center
        self.constraint = constraint
        self.repulsion_strength = repulsion_strength

        self.pos3d = {}

    #Subgraphs

    def make_subgraph(self, whole_graph: nx.Graph):
        """
        Subgraph from big PPI graph for this location with nodes with certain minimum of degrees (min_degree)
        """
        
        nodes_here = [n for n in self.node_ids if n in whole_graph]
        sg = whole_graph.subgraph(nodes_here).copy()
        sg.remove_nodes_from([n for n, d in sg.degree() if d <= self.min_degree]) # <- adjust degree here
        return sg

    # Simple metrics on Subgraphs

    def num_nodes(self, G: nx.Graph):
        H = self.make_subgraph(G)
        return H.number_of_nodes()

    def num_edges(self, G: nx.Graph):
        H = self.make_subgraph(G)
        return H.number_of_edges()
    

    def assign_positions_to_nodes(self, whole_graph: nx.Graph):
        """
        -Assigns all nodes of this location a (x,y,z) position (pos3D)
        -Saves them as Node Attributes in big whole-cell graph.
        
        """

        H = self.make_subgraph(whole_graph)
        
        self.pos3d = spring_layout_3d_constrained(H, self.constraint, seed=7, iterations=350, center=self.center, repulsion_strength = self.repulsion_strength)
        print(self.name, self.center)

        #set Node attribute in big Graph
        for node in self.pos3d.keys():
            whole_graph.nodes[node]["locations"][self.name] = self.pos3d[node]
        

class CellModel:
    """
    holds the whole_cell_graph and a ditionary of all location objects.
    Everything should come together here.
    """

    def __init__(self, whole_cell_graph: nx.Graph):
        self.whole_cell_graph = whole_cell_graph
        self.locations: Dict[str, Location] = {}

    def add_location(self, location: Location):
        self.locations[location.name] = location

    def assign_all_positions(self):
        """
        calls the placement function for all locations
        """
        for loc in self.locations.values():
            loc.assign_positions_to_nodes(self.whole_cell_graph)

    def summary(self):
        print(f"Whole-cell graph: {self.whole_cell_graph.number_of_nodes()} nodes, "
              f"{self.whole_cell_graph.number_of_edges()} edges")
        for name, loc in self.locations.items():
            H = loc.make_subgraph(self.whole_cell_graph)
            print(f"- {name}: {H.number_of_nodes()} nodes, {H.number_of_edges()} edges")
    
    def build_global_pos3d(self):
        """
        Merge all Location.pos3d into one dict: node_id -> (x, y, z)

        Decide which Proteins should be kept. If a node is in multiple locations, keep the first position
        encountered -> need to change that
        """
        pos3d = {}
        for loc in self.locations.values():
            for node_id, coords in loc.pos3d.items():
                if node_id not in pos3d: #must be changed for multilocalized proteins    !!!!!!!!!!!!!!!!!!!!!!!
                    pos3d[node_id] = coords
        return pos3d

    def plot_all_locations_3d(self, title="3D cell - all locations", show_boundaries=True):
        fig = go.Figure()

        def _translate_trace(trace, center):
            dx, dy, dz = center

            # Surface (2D-Gitter)
            if getattr(trace, "type", None) == "surface":
                trace.x = (np.asarray(trace.x) + dx)
                trace.y = (np.asarray(trace.y) + dy)
                trace.z = (np.asarray(trace.z) + dz)
                return trace

            # Scatter3d (1D-Listen, ggf. mit None als Trenner)
            def shift_1d(arr, d):
                return [None if v is None else v + d for v in arr]

            trace.x = shift_1d(trace.x, dx)
            trace.y = shift_1d(trace.y, dy)
            trace.z = shift_1d(trace.z, dz)
            return trace

        for loc_name, loc in self.locations.items():
            H = loc.make_subgraph(self.whole_cell_graph)

            if show_boundaries and hasattr(loc.constraint, "boundary_traces"):
                for tr in loc.constraint.boundary_traces():
                    tr = _translate_trace(tr, loc.center)
                    tr.name = f"{loc_name}: {getattr(tr, 'name', 'boundary')}"
                    tr.showlegend = False
                    fig.add_trace(tr)

            # Edges
            if H.number_of_edges() > 0:
                edge_x, edge_y, edge_z = [], [], []
                for u, v in H.edges():
                    x0, y0, z0 = loc.pos3d[u]
                    x1, y1, z1 = loc.pos3d[v]
                    edge_x += [x0, x1, None]
                    edge_y += [y0, y1, None]
                    edge_z += [z0, z1, None]

                fig.add_trace(go.Scatter3d(
                    x=edge_x, y=edge_y, z=edge_z,
                    mode="lines", line=dict(width=1),
                    hoverinfo="none", opacity=0.35,
                    showlegend=False
                ))

            # Nodes
            xs, ys, zs, texts = [], [], [], []
            for node_id, (x, y, z) in loc.pos3d.items():
                xs.append(x); ys.append(y); zs.append(z)
                texts.append(f"{node_id} ({loc_name})")

            fig.add_trace(go.Scatter3d(
                x=xs, y=ys, z=zs,
                mode="markers", name=loc_name,
                hovertext=texts, hoverinfo="text",
                marker=dict(size=2, opacity=0.35)
            ))

        fig.update_layout(
            title=title,
            showlegend=True,
            margin=dict(l=0, r=0, b=0, t=40),
            scene=dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False),
                aspectmode="data"
            ),
        )
        fig.show()



from whole_cell_spatial_graph import whole_cell_G, loc_to_proteins 
# -> where should the Graph come from?
# !!!!!!!!

if __name__ == "__main__":

    cell = CellModel(whole_cell_G)

    #Define Locations
    
    cytosol = Location(
        name="Cytosol",
        node_ids=loc_to_proteins["Cytosol"],
        min_degree=10,
        center=(0.0, 0.0, 0.0),
        constraint=EllipsoidConstraint(axes=(10000.0, 8000.0, 8000.0), wall=5000)
    )
    

    endoplasmicreticulum = Location(
        name="Endoplasmic reticulum",
        node_ids=loc_to_proteins["Endoplasmic reticulum"],
        min_degree=10,
        center=(0.0, 0.0, 0.0),
        constraint=ShellConstraint(inner_radius=4040.0, outer_radius=5000.0, wall=5000)
    )

    golgiapparatus = Location(
        name="Golgi apparatus",
        node_ids=loc_to_proteins["Golgi apparatus"],
        min_degree=10,
        center=(0.0, 0.0, 6000.0),
        constraint=EllipsoidConstraint(axes=(1500.0, 600.0, 400.0), wall=5000)
    )

    centrosome = Location(
        name="Centrosome",
        node_ids=loc_to_proteins["Centrosome"],
        min_degree=10,
        center=(0.0, 0.0, 5200.0),
        constraint=CylinderConstraint(radius=115.0, height=500.0, wall=5000)
    )

    intermediatefilaments = Location(
        name="Intermediate filaments",
        node_ids=loc_to_proteins["Intermediate filaments"],
        min_degree=0,
        center=(0.0, 0.0, 0.0),
        constraint=EllipsoidConstraint(axes=(9500.0, 7500.0, 7500.0), wall=5000)
    )

    microtubules = Location(
        name="Microtubules",
        node_ids=loc_to_proteins["Microtubules"],
        min_degree=0,
        center=(0.0, 0.0, 0.0),
        constraint=EllipsoidConstraint(axes=(9500.0, 7500.0, 7500.0), wall=5000)
    )

    actinfilaments = Location(
        name="Actin filaments",
        node_ids=loc_to_proteins["Actin filaments"],
        min_degree=10,
        center=(0.0, 0.0, 0.0),
        constraint=EllipsoidShellConstraint(axes=(9800.0, 7800.0, 7800.0), outer=1.020408, wall=500000)
    )

    plasmamembrane = Location(
        name="Plasma membrane",
        node_ids=loc_to_proteins["Plasma membrane"],
        min_degree=10,
        center=(0.0, 0.0, 0.0),
        constraint=EllipsoidShellConstraint(axes=(10000.0, 8000.0, 8000.0), outer=1.1, wall=500000)
    )

    nucleoplasm = Location(
        name="Nucleoplasm",
        node_ids=loc_to_proteins["Nucleoplasm"],
        min_degree=50,
        center=(0.0, 0.0, 0.0),
        constraint=SphereConstraint(radius=4000.0, wall=5000),
        repulsion_strength = 7
    )

    nuclearmembrane = Location(
        name="Nuclear membrane",
        node_ids=loc_to_proteins["Nuclear membrane"],
        min_degree=0,
        center=(0.0, 0.0, 0.0),
        constraint=ShellConstraint(inner_radius=4000.0, outer_radius=4040.0, wall=5000)
    )

    mitochondria = Location(
        name="Mitochondria",
        node_ids=loc_to_proteins["Mitochondria"],
        min_degree=10,
        center=(-5500.0, 2000.0, 0.0),
        constraint=EllipsoidConstraint(axes=(325.0, 325.0, 650.0), wall=5000)
    )

    primarycilium = Location(
        name="Primary cilium",
        node_ids=loc_to_proteins["Primary cilium"],
        min_degree=10,
        center=(0, 0, 8500),
        constraint=CylinderConstraint(radius=125.0, height=1000.0, wall=5000)
    )

    nucleoli = Location(
        name="Nucleoli",
        node_ids=loc_to_proteins["Nucleoli"],
        min_degree=10,
        center=(900.0, 800.0, 200.0),
        constraint=SphereConstraint(radius=600.0, wall=5000),
        repulsion_strength = 2
    )



    #Add locations to cell
    
    #cell.add_location(cytosol)
    
    cell.add_location(centrosome)
    cell.add_location(actinfilaments)
    cell.add_location(endoplasmicreticulum)
    cell.add_location(golgiapparatus)
    cell.add_location(intermediatefilaments)
    cell.add_location(microtubules)
    cell.add_location(nuclearmembrane)
    cell.add_location(nucleoli)
    cell.add_location(plasmamembrane)
    cell.add_location(nucleoplasm)
    cell.add_location(primarycilium)
    cell.add_location(mitochondria)
 
    print(cell.summary()) #prints a summary of all locations

    cell.assign_all_positions() #assigns 3D coordinates to locations

    #print(primarycilium.pos3d)

    '''
    from constraint_layout_classes import *

    test_subgraph = mitochondria.make_subgraph(whole_cell_G)

    constraint = SphereConstraint(radius=5.0, wall=5000)
    pos = spring_layout_3d_constrained(test_subgraph, constraint, seed=7, iterations=450)
    print(pos)
    plot_layout_3d(test_subgraph, pos, extra_traces=constraint.boundary_traces(), title="Mitochondria").show()
    '''

    #Print Edge list of whole_cell_graph in cell
    '''
    edge_nodes_with_attrs = [
        ((u, dict(cell.whole_cell_graph.nodes[u])), (v, dict(cell.whole_cell_graph.nodes[v])))
        for u, v in cell.whole_cell_graph.edges()
    ]

    for pair in edge_nodes_with_attrs[:10]:
        #print(pair)
        continue
    '''
    #Plot 
    cell.plot_all_locations_3d("3D cell - all locations with edges")
