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
    - center:           3D coordinates of center (x, y, z)
    - radius:           Size parameter
    - sample_function:  3D coordinate sampling geometry
    """

    def __init__(
        self,
        name: str,
        node_ids: List[str],
        center: Tuple[float, float, float],
        radius: float,
        geometry: Callable,
        constraint: object
    ):
        self.name = name
        self.node_ids = node_ids
        self.center = center
        self.radius = radius
        self.pos3d = {}
        self.geometry = geometry
        self.constraint = constraint

    #Subgraphs

    def make_subgraph(self, whole_graph: nx.Graph):
        """Subgraph from big PPI graph for this location with a certain degree."""
        nodes_here = [n for n in self.node_ids if n in whole_graph]
        sg = whole_graph.subgraph(nodes_here).copy()
        sg.remove_nodes_from([n for n, d in sg.degree() if d < 1]) # <- adjust degree here
        return sg

    # Simple metrics on Subgraphs

    def num_nodes(self, G: nx.Graph):
        H = self.make_subgraph(G)
        return H.number_of_nodes()

    def num_edges(self, G: nx.Graph):
        H = self.make_subgraph(G)
        return H.number_of_edges()
    
    """
    def assign_positions_to_nodes(self, whole_graph: nx.Graph):
        
        #-Assigns all nodes of this location a (x,y,z) position (pos3D)
        #-Saves them as Node Attributes in big whole-cell graph.
        
        H = self.make_subgraph(whole_graph)

        layout = nx.spring_layout(H, dim=3, seed=42) #creates a force-directed layout: node_id -> np.array([x, y, z])
        #print(layout)
        
        coords = np.array(list(layout.values())) #get only coordinates
        center_layout = coords.mean(axis=0) #calculate middle point of coordinates

        for node in H.nodes():
            # Direction from Layout
            v = np.array(layout[node]) - center_layout #vector from network center to this node
            norm = np.linalg.norm(v) #length of vector
            
            if norm < 1e-9: #if vector too small, because point is too similar to middle-point: use fallback
                v = np.array([1.0, 0.0, 0.0])  #simple vector (change?)

            direction = v / np.linalg.norm(v) #Calculate Einheitsvector (lenght of 1): gives direction

            #print(np.linalg.norm(direction))

            #get point on surface by projecting direction
            x, y, z = self.geometry(direction, self.center, self.radius)

            #set Node attribute in big Graph
            whole_graph.nodes[node]["locations"][self.name] = [x, y, z]

            #save node coordinate in location positions
            self.pos3d[node] = [x, y, z]
        
        print(self.pos3d)
    """

    def assign_positions_to_nodes(self, whole_graph: nx.Graph):
        """
        -Assigns all nodes of this location a (x,y,z) position (pos3D)
        -Saves them as Node Attributes in big whole-cell graph.
        
        """

        H = self.make_subgraph(whole_graph)
        
        self.pos3d = spring_layout_3d_constrained(H, self.constraint, seed=7, iterations=450, center=self.center)
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
                if node_id not in pos3d: #must be changed for multilocalized proteins
                    pos3d[node_id] = coords
        return pos3d

    def plot_all_locations_3d(self, title = "3D cell - all locations"):
        fig = go.Figure()

        
        for loc_name, loc in self.locations.items():    # Loop over all Locations
            
            # Subgraph with only nodes of one location
            H = loc.make_subgraph(self.whole_cell_graph)

            # create Edges within a location
            if H.number_of_edges() > 0:
                edge_x, edge_y, edge_z = [], [], []

                for u, v in H.edges():

                    x0, y0, z0 = loc.pos3d[u]
                    x1, y1, z1 = loc.pos3d[v]

                    edge_x += [x0, x1, None] #Plotly uses None as separation mark
                    edge_y += [y0, y1, None]
                    edge_z += [z0, z1, None]

                fig.add_trace(
                    go.Scatter3d(
                        x=edge_x,
                        y=edge_y,
                        z=edge_z,
                        mode="lines",
                        line=dict(width=1),
                        hoverinfo="none",
                        opacity=0.35,
                        showlegend=False,  # Edges not in Legend
                    )
                )

            # nodes of this location

            xs, ys, zs, texts = [], [], [], []
            for node_id, (x, y, z) in loc.pos3d.items():
                xs.append(x)
                ys.append(y)
                zs.append(z)
                texts.append(f"{node_id} ({loc_name})")

            fig.add_trace(
                go.Scatter3d(
                    x=xs,
                    y=ys,
                    z=zs,
                    mode="markers",
                    name=loc_name,      #Nodes in legend
                    hovertext=texts,
                    hoverinfo="text",
                    marker=dict(        #Visual properties
                        size=2,
                        opacity=0.35,
                    ),
                )
            )

        fig.update_layout(
            title=title,
            showlegend=True,
            margin=dict(l=0, r=0, b=0, t=40),
            scene=dict(
                xaxis=dict(visible=False),  #Hide all UI axes
                yaxis=dict(visible=False),
                zaxis=dict(visible=False),
                aspectmode="data" 
            ),
        )

        fig.show()


from whole_cell_spatial_graph import whole_cell_G, loc_to_proteins 
# -> where should the Graph come from?

cell = CellModel(whole_cell_G)

#Define Locations
'''
cytosol = Location(
    name="Cytosol",
    node_ids=loc_to_proteins["Cytosol"],
    center=(50.0, 0.0, 0.0),
    radius=(9, 2, 1),
    geometry=project_to_sphere
)

centrosome = Location(
    name="Centrosome",
    node_ids=loc_to_proteins["Centrosome"],
    center=(100.0, 0.0, 0.0),
    radius=4.0,
    geometry=project_to_sphere
)

endoplasmicreticulum = Location(
    name="Endoplasmic reticulum",
    node_ids=loc_to_proteins["Endoplasmic reticulum"],
    center=(0.0, 50.0, 0.0),
    radius=20.0,
    geometry=project_to_sphere
)

golgiapparatus = Location(
    name="Golgi apparatus",
    node_ids=loc_to_proteins["Golgi apparatus"],
    center=(50.0, 50.0, 0.0),
    radius=10.0,
    geometry=project_to_sphere
)

intermediatefilaments = Location(
    name="intermediatefilaments",
    node_ids=loc_to_proteins["Intermediate filaments"],
    center=(50.0, 50.0, 50.0),
    radius=4.0,
    geometry=project_to_sphere
)

microtubules = Location(
    name="Microtubules",
    node_ids=loc_to_proteins["Microtubules"],
    center=(50.0, 50.0, 100.0),
    radius=(50, 50, 50),
    geometry=project_to_sphere
)

nuclearmembrane = Location(
    name="Nuclear membrane",
    node_ids=loc_to_proteins["Nuclear membrane"],
    center=(150.0, 150.0, 150.0),
    radius=1.0,
    geometry=project_to_sphere
)

plasmamembrane = Location(
    name="Plasma membrane",
    node_ids=loc_to_proteins["Plasma membrane"],
    center=(0.0, 0.0, 0.0),
    radius=900.0,
    geometry=project_to_sphere_surface
)


actinfilaments = Location(
    name="Actin filaments",
    node_ids=loc_to_proteins["Actin filaments"],
    center=(0.0, 0.0, 0.0),
    radius=500.0,
    geometry=sample_sphere_surface_projection
)

nucleoplasm = Location(
    name="Nucleoplasm",
    node_ids=loc_to_proteins["Nucleoplasm"],
    center=(0, 0, 0),
    radius=9.0,
    geometry=sample_ellipsoid_volume_projection,
    constraint=SphereConstraint(radius=5.0, wall=5000)

)

'''
mitochondria = Location(
    name="Mitochondria",
    node_ids=loc_to_proteins["Mitochondria"],
    center=(0.0, 20.0, 0.0),
    radius=(0, 100, 0),
    geometry=sample_ellipsoid_volume_projection,
    constraint=EllipsoidConstraint(axes=(3.25, 3.25, 6.50), wall=5000)
)

primarycilium = Location(
    name="Primary cilium",
    node_ids=loc_to_proteins["Primary cilium"],
    center=(0.0, 0.0, 0.0),
    radius=250.0,
    geometry=sample_sphere_surface_projection,
    constraint=ShellConstraint(inner_radius=5.1, outer_radius=5.3, wall=5000)
)


nucleoli = Location(
    name="Nucleoli",
    node_ids=loc_to_proteins["Nucleoli"],
    center=(0.0, 0.0, 0.0),
    radius=9.0,
    geometry=sample_ellipsoid_volume_projection,
    constraint=SphereConstraint(radius=5.0, wall=5000)
)

#Add locations to cell
'''
cell.add_location(cytosol)
cell.add_location(centrosome)
cell.add_location(actinfilaments)
cell.add_location(endoplasmicreticulum)
cell.add_location(golgiapparatus)
cell.add_location(intermediatefilaments)
cell.add_location(microtubules)
cell.add_location(mitochondria)
cell.add_location(nuclearmembrane)
cell.add_location(nucleoli)
cell.add_location(actinfilaments)
cell.add_location(primarycilium)
cell.add_location(nucleoli)
cell.add_location(nucleoplasm)
'''
cell.add_location(primarycilium)
#cell.add_location(mitochondria)
cell.add_location(nucleoli)

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
