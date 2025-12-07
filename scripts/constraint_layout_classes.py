
import numpy as np
import networkx as nx
import math
import plotly.graph_objects as go
from typing import Optional, List

"""
Constraint Class:
    - sample()          ... sample initial start positions for layout
    - forces()          ... calculate forces in right direction with appropiate strength for each point q
    - boundary_traces() ... draw boundaries as 3D plotly GO surface
"""
class SphereConstraint():
    def __init__(self, radius: float = 1.0, wall: float = 5000.0):
        self.radius = float(radius)
        self.wall = float(wall)

    def sample(self, n, rng):
        directions = rng.normal(size=(n, 3))  # Generates random direction vectors with a (0,0,0) center
        directions /= (np.linalg.norm(directions, axis=1, keepdims=True) + 1e-12) # normalizes all these vectors, so they all have a length of 1 + protection against division by 0
        length = rng.random(n) ** (1.0 / 3.0) #randomly select lengths for vectors between 0 and 1 -> volume. 1/3 for equal distribution in 3D  

        return directions * (length * self.radius)[:, None]
        #(n,) -> [:, None] -> (n,1) 
        #-> multiply radius with lengths and directions

    def forces(self, q):
        """
        q           ... all current point positions
        returns F   ... force vector for each point
        """
        radius = self.radius

        dist = np.linalg.norm(q, axis=1) + 1e-12 #scalar: Distance of each point q to center

        directions = q / dist[:, None] 
            #direction of points as unit vectors
            #length of 1 important to build radial forces -> separate force strength an direction
            #alway points "outside"

        # soft wall (outside only)
        excess = np.maximum(0.0, dist - radius)
            # if point inside of sphere: excess = 0
            # if point outside of sphere:  excess = dist - radius

        F = -(self.wall * excess)[:, None] * directions # calculate force for each point q
            # minus for direction inside
            # wall is an amplification factor for force strength in right direction
            # only active when point outside of sphere

        return F

    #3D Boundary Shape (for Plotly)
    def boundary_traces(self):
        radius = self.radius
        uu = np.linspace(0, 2*np.pi, 45)
        vv = np.linspace(0, np.pi, 24)
        xs = radius * np.outer(np.cos(uu), np.sin(vv))
        ys = radius * np.outer(np.sin(uu), np.sin(vv))
        zs = radius * np.outer(np.ones_like(uu), np.cos(vv))
        return [go.Surface(x=xs, y=ys, z=zs, opacity=0.08, showscale=False,
                           hoverinfo="skip", name="boundary")]


class ShellConstraint():
    # Spherical shell with inner/outer radius.
    def __init__(self, inner_radius: float = 0.6, outer_radius: float = 1.0, wall: float = 5000.0):
        self.inner_radius = float(inner_radius)
        self.outer_radius = float(outer_radius)
        self.wall = float(wall)

    def sample(self, n, rng):
        #create random vectors and norm then -> random directions on unit-sphere (r=1)
        directions = rng.normal(size=(n, 3))
        directions /= (np.linalg.norm(directions, axis=1, keepdims=True) + 1e-12)
        
        #set length to be in intervall between inner and outer radius and shift everything by inner radius
        length = (self.inner_radius + rng.random(n) * (self.outer_radius - self.inner_radius))
        return directions * length[:, None]
            #(n,) -> [:, None] -> (n,1) 

    def forces(self, q):
        rin = self.inner_radius
        rout = self.outer_radius
        
        dist = np.linalg.norm(q, axis=1) + 1e-12 #scalar: Distance of each point q to center

        directions = q / dist[:, None]
            #direction of points as unit vectors
            #length of 1 important to build radial forces -> separate force strength an direction
            #alway points "outside"

        excess_out = np.maximum(0.0, dist - rout) #Force Strength in
        excess_in = np.maximum(0.0, rin - dist) #Force Strength out

        F = np.zeros_like(q)
        F += -(self.wall * excess_out)[:, None] * directions  # pull in
        F += +(self.wall * excess_in)[:, None] * directions   # push out

        return F

    #3D Boundary Shape (for Plotly)
    def boundary_traces(self):
        def sphere_surface(R, opacity, name):
            uu = np.linspace(0, 2*np.pi, 45)
            vv = np.linspace(0, np.pi, 24)
            xs = R * np.outer(np.cos(uu), np.sin(vv))
            ys = R * np.outer(np.sin(uu), np.sin(vv))
            zs = R * np.outer(np.ones_like(uu), np.cos(vv))
            return go.Surface(x=xs, y=ys, z=zs, opacity=opacity, showscale=False,
                              hoverinfo="skip", name=name)

        return [
            sphere_surface(self.outer_radius, 0.06, "outer"),
            sphere_surface(self.inner_radius, 0.10, "inner"),
        ]


class CylinderConstraint():
    def __init__(self, radius: float = 1.0, height: float = 2.0, wall: float = 5000.0):
        self.radius = float(radius)
        self.height = float(height)
        self.wall = float(wall)

    def sample(self, n, rng):
        #get random directions on unit circle in xy-plane
        theta = rng.random(n) * 2 * np.pi #(n, )
        #print(theta)

        #create random radial 2D directions via theta -> only x and y
        directions_xy = np.column_stack([np.cos(theta), np.sin(theta)]) 
            # np.cos/sin(theta) creates 1D (n,) array
            # column_stack glues these 1D array into columns -> (n,2)
        
        radial_length = np.full(n, self.radius)  #(n,) -> radial length array, where all entries = self.radius

        xy = directions_xy * radial_length[:, None]  #(n,2)
            #(n,) -> [:, None] -> (n,1) 

        #get random z via height
        z = rng.random(n) * self.height
        
        return np.column_stack([xy, z]) #glue x,y and z togehter in one (n,3) array

    def forces(self, q):
        radius = self.radius
        half = 0.5 * self.height #distance from center to caps (z-coordinate)

        #unpack each coordinate
        x, y, z = q[:, 0], q[:, 1], q[:, 2]

        dist_radial = np.sqrt(x*x + y*y) + 1e-12
            #shortest perpendicular radial distance to the z-axis + divison-protection
            # length of a vector is sqr(x²+y²+z²) -> ignore z to know how far point is away from z-axis = radial distance

        #get normalized direction (length = 1 in xy plane) which points away from z-axis in the direction of the mantle of the cylinder
        # only in xy-plane -> z = zero
        direction_radial = np.column_stack([x/dist_radial, y/dist_radial, np.zeros_like(dist_radial)])

        #wall forces
        # outside mantle
        excess_radial = np.maximum(0.0, dist_radial - radius)
        F = -(self.wall * excess_radial)[:, None] * direction_radial

        # outside caps
        excess_z = np.maximum(0.0, np.abs(z) - half)         
        F[:, 2] += -(self.wall * excess_z) * np.sign(z)
            #F[:, 2] -> set force of z-coordinates
            #determine direction of cap force via np.sign -> push up or down

        return F

    def boundary_traces(self):
        #build 3 surfaces: Top cap, bottom cap and mantle
        R, H = self.radius, self.height
        half = H / 2.0

        theta = np.linspace(0, 2*np.pi, 100)
        cx = R * np.cos(theta)
        cy = R * np.sin(theta)

        top = go.Scatter3d(x=cx, y=cy, z=np.full_like(cx, half), mode="lines",
                           line=dict(width=2, color="rgba(160,160,160,0.35)"),
                           hoverinfo="skip", name="boundary")
        bot = go.Scatter3d(x=cx, y=cy, z=np.full_like(cx, -half), mode="lines",
                           line=dict(width=2, color="rgba(160,160,160,0.35)"),
                           hoverinfo="skip", showlegend=False)

        vl_x, vl_y, vl_z = [], [], []
        for th in np.linspace(0, 2*np.pi, 12, endpoint=False):
            vl_x += [R*np.cos(th), R*np.cos(th), None]
            vl_y += [R*np.sin(th), R*np.sin(th), None]
            vl_z += [-half, half, None]
        side = go.Scatter3d(x=vl_x, y=vl_y, z=vl_z, mode="lines",
                            line=dict(width=1, color="rgba(160,160,160,0.25)"),
                            hoverinfo="skip", showlegend=False)

        return [top, bot, side]

class EllipsoidConstraint():
    """
    Solid Ellipsoid: (x/a)^2 + (y/b)^2 + (z/c)^2 <= 1

    - sample(): uniform in ellipsoid volume via "unit ball sample" -> scale by axes
    - forces(): soft wall in ellipsoid metric (only outside)
    - boundary_traces(): plotly surface for the ellipsoid boundary
    """
    def __init__(self, axes=(3.0, 5.0, 3.0), wall: float = 5000.0):
        if len(axes) != 3:
            raise ValueError("axes must be a 3-tuple (a,b,c)")
        a, b, c = float(axes[0]), float(axes[1]), float(axes[2])
        if a <= 0 or b <= 0 or c <= 0:
            raise ValueError("axes values must be > 0")
        self.axes = (a, b, c)
        self.wall = float(wall)
        self.scale = self.axes

    def sample(self, n, rng):
        # uniform in unit ball
        directions = rng.normal(size=(n, 3))
        directions /= (np.linalg.norm(directions, axis=1, keepdims=True) + 1e-12)
        length = rng.random(n) ** (1.0 / 3.0)
        p = directions * length[:, None]  # in unit ball

        # scale into ellipsoid volume
        a, b, c = self.axes
        p[:, 0] *= a
        p[:, 1] *= b
        p[:, 2] *= c
        return p

    def forces(self, q):
        a, b, c = self.axes
        axes = np.array([a, b, c], dtype=float)
        axes2 = axes * axes

        # ellipsoid "distance" in scaled space: s = ||(x/a, y/b, z/c)||
        p = q / axes[None, :]
        s = np.linalg.norm(p, axis=1) + 1e-12

        # direction for pushing back inward in ORIGINAL coordinates:
        # grad(s) w.r.t q is (q/axes^2)/s
        directions = (q / axes2[None, :]) / s[:, None]

        # soft wall only outside (s>1)
        excess = np.maximum(0.0, s - 1.0)
        F = -(self.wall * excess)[:, None] * directions
        return F

    def boundary_traces(self):
        a, b, c = self.axes
        uu = np.linspace(0, 2*np.pi, 45)
        vv = np.linspace(0, np.pi, 24)
        xs = a * np.outer(np.cos(uu), np.sin(vv))
        ys = b * np.outer(np.sin(uu), np.sin(vv))
        zs = c * np.outer(np.ones_like(uu), np.cos(vv))
        return [go.Surface(x=xs, y=ys, z=zs, opacity=0.08, showscale=False,
                           hoverinfo="skip", name="boundary")]

class EllipsoidShellConstraint():
    """
    Ellipsoid shell: 1 <= s <= outer, where s = ||(x/a, y/b, z/c)|| in scaled space.
    Equivalently: inner surface is the ellipsoid with axes (a,b,c),
                  outer surface is the ellipsoid with axes (outer*a, outer*b, outer*c).

    - sample(): uniform in shell volume via unit-ball shell sample -> scale by axes
    - forces(): soft walls at inner (s<1) and outer (s>outer)
    """
    def __init__(self, axes=(3.0, 5.0, 3.0), outer: float = 1.3, wall: float = 5000.0):
        if len(axes) != 3:
            raise ValueError("axes must be a 3-tuple (a,b,c)")
        a, b, c = float(axes[0]), float(axes[1]), float(axes[2])
        if a <= 0 or b <= 0 or c <= 0:
            raise ValueError("axes values must be > 0")
        outer = float(outer)
        if outer <= 1.0:
            raise ValueError("outer must be > 1.0 (shell thickness)")
        self.axes = (a, b, c)
        self.outer = outer
        self.wall = float(wall)
        # useful for anything that wants a rough scale
        self.scale = (outer*a, outer*b, outer*c)

    def sample(self, n, rng):
        # sample uniform in unit BALL SHELL: s in [1, outer] with volume-uniform radius
        directions = rng.normal(size=(n, 3))
        directions /= (np.linalg.norm(directions, axis=1, keepdims=True) + 1e-12)

        inner3 = 1.0 ** 3
        outer3 = self.outer ** 3
        s = (inner3 + rng.random(n) * (outer3 - inner3)) ** (1.0 / 3.0)  # volume-uniform shell radius

        p = directions * s[:, None]  # points in unit shell (scaled space)

        # scale into ellipsoid shell
        a, b, c = self.axes
        p[:, 0] *= a
        p[:, 1] *= b
        p[:, 2] *= c
        return p

    def forces(self, q):
        a, b, c = self.axes
        axes = np.array([a, b, c], dtype=float)
        axes2 = axes * axes

        # scaled radius s = ||(x/a, y/b, z/c)||
        p = q / axes[None, :]
        s = np.linalg.norm(p, axis=1) + 1e-12

        # gradient direction in original coordinates: grad(s) = (q/axes^2)/s
        directions = (q / axes2[None, :]) / s[:, None]

        # outer wall: s > outer -> push inward
        excess_out = np.maximum(0.0, s - self.outer)
        # inner wall: s < 1 -> push outward
        excess_in = np.maximum(0.0, 1.0 - s)

        F = np.zeros_like(q)
        F += -(self.wall * excess_out)[:, None] * directions
        F += +(self.wall * excess_in)[:, None] * directions
        return F

    def boundary_traces(self):
        a, b, c = self.axes
        uu = np.linspace(0, 2*np.pi, 45)
        vv = np.linspace(0, np.pi, 24)

        def ell_surface(ax, by, cz, opacity, name):
            xs = ax * np.outer(np.cos(uu), np.sin(vv))
            ys = by * np.outer(np.sin(uu), np.sin(vv))
            zs = cz * np.outer(np.ones_like(uu), np.cos(vv))
            return go.Surface(x=xs, y=ys, z=zs, opacity=opacity, showscale=False,
                              hoverinfo="skip", name=name)

        outer = self.outer
        return [
            ell_surface(outer*a, outer*b, outer*c, 0.06, "outer"),
            ell_surface(a, b, c, 0.10, "inner"),
        ]


def spring_layout_3d_constrained(
    G,
    constraint,
    iterations: int = 350,
    seed: int = 1,
    k: Optional[float] = None,
    threshold: float = 1e-4,
    recenter_each_iter: bool = True,
    center=(0.0, 0.0, 0.0),
    edge_strength: float = 1.0,      # <--- Attraktion
    repulsion_strength: float = 1.0, # <--- Repulsion 
):
    nodes = list(G)
    n = len(nodes)

    #Adjazenzmatrix: 1 if edge, 0 if no edge
    A = nx.to_numpy_array(G, nodelist=nodes, weight=None, dtype=float)

    rng = np.random.default_rng(seed) #reprodcueable random generator
    pos = constraint.sample(n, rng) #sample starting positions

    # k automatisch an Container-Skala koppeln
    if k is None:
        extent = np.ptp(pos, axis=0)              # (Lx, Ly, Lz)
        extent = np.maximum(extent, 1e-6)         # avoid zeros
        bbox_vol = float(extent[0] * extent[1] * extent[2])
        k = (bbox_vol / n) ** (1.0 / 3.0)         # ~typical spacing in 3D

        k *= 0.8
        print(k)

    #Temperature
    t = np.ptp(pos, axis=0).max() * 0.1
    t = 0.1 if t == 0 else t #at least 0.1
    dt = t / (iterations + 1) #reduction of temperature per interation

    for _ in range(iterations):
        #Create Matrix of all point combinations (n_i, n_j, 3) with direction vectors from i to j
        delta = pos[:, None, :] - pos[None, :, :]
        
        #Create Matrix of all point combinations (n_i, n_j, 3) with scalar distances from i to j
        #||pos[i] - pos[j]||
        dist = np.linalg.norm(delta, axis=-1)
        
        np.clip(dist, 0.01, None, out=dist) 
            #avoid dist to be smaller then 0.01 -> no 0 divison
            #e.g. when i = j -> dist is 0

        # Fruchterman-Reingold coefficients -> repulsion minus attraction
        # repulsion ~ k^2/dist^2
        # attraction ~ dist/k (only edges)
        
        # coeff is a n x n matrix of scalars, which determines for each point how strong they attract or repell
        coeff = (repulsion_strength * (k * k) / (dist * dist)) - (edge_strength * A * dist / k)
        
        # disp is for each node i the whole move-vector (force), which is calculated from all other nodes
        disp = np.einsum("ijk,ij->ik", delta, coeff)

        disp += constraint.forces(pos) # Constraint add force

        length = np.linalg.norm(disp, axis=1)
        length = np.clip(length, 0.01, None)
        step = disp * (t / length)[:, None]
        pos += step

        if recenter_each_iter:
            pos -= pos.mean(axis=0)

        t -= dt
        if (np.linalg.norm(step) / n) < threshold:
            break

    #print(type(pos), getattr(pos, "shape", None))
    #print(type(center), np.asarray(center).shape)   
    
    pos += np.asarray(center, dtype=float)
    
    return {nodes[i]: pos[i] for i in range(n)}


"""
def plot_layout_3d(G, pos, extra_traces=None, title="3D constrained spring layout"):
    #Plot only nodes + edges + boundaries (as extra_trace)
    extra_traces = extra_traces or []
    nodes = list(G)

    X = np.array([pos[u][0] for u in nodes])
    Y = np.array([pos[u][1] for u in nodes])
    Z = np.array([pos[u][2] for u in nodes])

    xe, ye, ze = [], [], []
    for u, v in G.edges():
        xe += [pos[u][0], pos[v][0], None]
        ye += [pos[u][1], pos[v][1], None]
        ze += [pos[u][2], pos[v][2], None]

    edge_trace = go.Scatter3d(
        x=xe, y=ye, z=ze, mode="lines",
        line=dict(width=2, color="rgba(120,120,120,0.35)"),
        hoverinfo="skip", name="edges",
    )

    C = np.column_stack([X, Y, Z]).mean(axis=0)
    cval = np.sqrt((X - C[0])**2 + (Y - C[1])**2 + (Z - C[2])**2)

    node_trace = go.Scatter3d(
        x=X, y=Y, z=Z, mode="markers",
        marker=dict(size=6, color=cval, colorscale="Viridis", opacity=0.95,
                    colorbar=dict(title="dist to centroid")),
        hoverinfo="text", name="nodes",
        text=[f"node {u}<br>deg={G.degree(u)}" for u in nodes],
    )

    fig = go.Figure(data=list(extra_traces) + [edge_trace, node_trace])
    
    fig.update_layout(
        title=title,
        showlegend=True,
        margin=dict(l=0, r=0, b=0, t=40),
        scene=dict(
            xaxis=dict(visible=False),  #Hide all UI axes
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            aspectmode="data",
        ),
    )
    return fig
"""


#G = nx.barabasi_albert_graph(800, 2, seed=2)

#sphere = SphereConstraint(radius=5.0, wall=5000)
#pos = spring_layout_3d_constrained(G, sphere, seed=7, iterations=400)
#plot_layout_3d(G, pos, extra_traces=sphere.boundary_traces(), title="Sphere").show()


#cyl = CylinderConstraint(radius=3, height=10, wall=5000)
#pos = spring_layout_3d_constrained(G, cyl, seed=7)
#plot_layout_3d(G, pos, extra_traces=cyl.boundary_traces(), title="Cylinder").show()

#shell = ShellConstraint(inner_radius=1, outer_radius=1.2, wall=5000)
#pos = spring_layout_3d_constrained(G, shell, seed=7, iterations=350)
#plot_layout_3d(G, pos, extra_traces=shell.boundary_traces(), title="Shell").show()

#ell = EllipsoidConstraint(axes=(6, 3, 10), wall=5000)
#pos = spring_layout_3d_constrained(G, ell, seed=7, iterations=400)
#plot_layout_3d(G, pos, extra_traces=ell.boundary_traces(), title="Ellipsoid").show()

#ell_shell = EllipsoidShellConstraint(axes=(10, 8, 8), outer=1.01, wall=500000)
#pos = spring_layout_3d_constrained(G, ell_shell, seed=7)
#plot_layout_3d(G, pos, extra_traces=ell_shell.boundary_traces(), title="Ellipsoid Shell").show()
