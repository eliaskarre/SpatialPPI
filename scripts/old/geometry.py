# src/geometry.py
# Minimal geometry generator stubs for organelle placement.
# Returns point clouds for simple parametric organelles.
import numpy as np
import math
import random


#random Point-Cloud sampling functions

def sample_sphere_cloud(center, radius, n_points=1000, rng=None):
    rng = np.random.default_rng(rng)
    u = rng.random(n_points)
    v = rng.random(n_points)
    theta = 2 * np.pi * u
    phi = np.arccos(2 * v - 1)
    x = radius * np.sin(phi) * np.cos(theta) + center[0]
    y = radius * np.sin(phi) * np.sin(theta) + center[1]
    z = radius * np.cos(phi) + center[2]
    return np.vstack((x, y, z)).T

def sample_ellipsoid_cloud(center, axes, n_points=1000, rng=None):
    #axes = (a, b, c)
    rng = np.random.default_rng(rng)
    u = rng.random(n_points)
    v = rng.random(n_points)
    theta = 2 * np.pi * u
    phi = np.arccos(2 * v - 1)
    a, b, c = axes
    x = a * np.sin(phi) * np.cos(theta) + center[0]
    y = b * np.sin(phi) * np.sin(theta) + center[1]
    z = c * np.cos(phi) + center[2]
    return np.vstack((x, y, z)).T

def sample_shell_surface_cloud(center, axes, thickness, n_points=1000, rng=None):
    # For membrane-proximal placements: sample two concentric ellipsoids and pick shell points.
    outer = sample_ellipsoid_cloud(center, axes, n_points, rng=rng)
    inner_axes = (max(axes[0] - thickness, 0.001), max(axes[1] - thickness, 0.001), max(axes[2] - thickness, 0.001))
    inner = sample_ellipsoid_cloud(center, inner_axes, n_points, rng=rng)
    return outer

#random single-point sampling functions

#Surface of sphere
def sample_sphere_point(center, radius):
    cx, cy, cz = center
    
    theta = 2 * np.pi * random.random()
    phi = np.arccos(2 * random.random() - 1)

    x = radius * np.sin(phi) * np.cos(theta) + cx
    y = radius * np.sin(phi) * np.sin(theta) + cy
    z = radius * np.cos(phi) + cz
    return x, y, z #one point

#Surface of ellipsoid
def sample_ellipsoid_point(center, axes):
    cx, cy, cz = center
    
    theta = 2 * np.pi * random.random()
    phi = np.arccos(2 * random.random() - 1)
    a, b, c = axes
    x = a * np.sin(phi) * np.cos(theta) + cx
    y = b * np.sin(phi) * np.sin(theta) + cy
    z = c * np.cos(phi) + cz
    
    return x, y, z

#Shell
def sample_shell_surface(center, axes, thickness=70):
    # For membrane-proximal placements: sample two concentric ellipsoids and pick shell points.
    outer = sample_ellipsoid_point(center, axes)
    inner_axes = (
        max(axes[0] - thickness, 0.001),
        max(axes[1] - thickness, 0.001),
        max(axes[2] - thickness, 0.001),
    )
    inner = sample_ellipsoid_point(center, inner_axes)
    return outer

# Sampling of projection from directions Functions

#Sphere Surface
def sample_sphere_surface_projection(direction, center, radius):
    d = direction / np.linalg.norm(direction) #normalize (input should be a einheitsvektorn, but just in case)
    cx, cy, cz = center
    x = cx + radius * d[0]
    y = cy + radius * d[1]
    z = cz + radius * d[2]
    return x, y, z

#Sphere Volume
def sample_sphere_volume_projection(direction, center, radius):
    d = direction / np.linalg.norm(direction) #normalize (input should be a einheitsvektorn, but just in case)
    cx, cy, cz = center

    # random radius in volume
    u = np.random.rand()
    r = radius * (u ** (1.0 / 3.0)) # 1/3 for u ~ Uniform(0,1) -> big radius have much more space then small ones
    x = cx + r * d[0]
    y = cy + r * d[1]
    z = cz + r * d[2]
    return x, y, z

#Cylinder Volume
def sample_cylinder_volume_projection(direction, center, radius, height = 60):
    """
    Use direction only to define the angular direction in the xy-plane.
    Sample radius & height such that the point is uniform in the cylinder volume.
    -> Cylinder axis = z-axis
    """
    cx, cy, cz = center
    d = np.asarray(direction, dtype=float)
    d = d / np.linalg.norm(d)

    # Project direction onto xy-plane
    d_xy = np.array([d[0], d[1]], dtype=float)  #only x and y of direction
    norm_xy = np.linalg.norm(d_xy)              #how far is direction in xy-plane

    # If direction is parallel to z-axis (almost no length in xy-plane) choose random angle
    if norm_xy < 1e-8:
        phi = 2 * np.pi * np.random.rand()
        d_xy = np.array([np.cos(phi), np.sin(phi)]) #sin and cos convert angle into x and y
    else:
        d_xy /= norm_xy     #normalize xy-direction to length of 1

    # sample xy-plane (radial distance)
    u_r = np.random.rand()
    r = radius * np.sqrt(u_r) #Square root for uniform sampling

    # get x and y from plane
    x = cx + r * d_xy[0]
    y = cy + r * d_xy[1]

    # get z (height)
    u_h = np.random.rand()
    z = cz + (u_h - 0.5) * height # -0.5 for uniform sampling

    return x, y, z

#Ellipsoid Volume
def sample_ellipsoid_volume_projection(direction, center, axes):
    """
    Sample a point inside the ellipsoid volume along a given direction.
    The direction is fixed; the distance from the center is random in [0, t_max]
    with r ∝ u^(1/3), analogous to uniform sampling in a ball (if directions
    are chosen uniformly overall).
    """

    #Convert to array an normalise
    d = np.asarray(direction, dtype=float)
    d = d / np.linalg.norm(d)

    #Unpack center and axes
    cx, cy, cz = center
    a, b, c = axes

    #distance along d to ellipsoid surface
    denom = (d[0]**2 / a**2) + (d[1]**2 / b**2) + (d[2]**2 / c**2)
    t_max = 1.0 / np.sqrt(denom) # # calculcates the length between center and ellipsiod-wall for a direction d

    # random point on length between center and ellipsiod-wall
    u = np.random.rand()
    r = t_max * (u ** (1.0 / 3.0)) # r^3 proportional to u, like in sphere for uniform sampling

    #get x, y, z coordinates
    x = cx + r * d[0]
    y = cy + r * d[1]
    z = cz + r * d[2]
    return x, y, z

#Shell Volume
'''
def sample_sphere_shell_volume_projection(direction, center, r_inner, r_outer = 700):
    """
    Sample a point uniformly in the spherical shell volume between r_inner and r_outer
    along a given direction from the center.
    direction: beliebiger Vektor (wird normalisiert)
    center: (cx, cy, cz)
    r_inner, r_outer: innerer und äußerer Radius der Schale
    """
    d = np.asarray(direction, dtype=float)
    d = d / np.linalg.norm(d)

    cx, cy, cz = center

    # Zufälliger Radius in Shell-Volumen (uniform im Volumen)
    u = np.random.rand()
    r_inner3 = r_inner ** 3
    r_outer3 = r_outer ** 3
    r = (r_inner3 + u * (r_outer3 - r_inner3)) ** (1.0 / 3.0)

    x = cx + r * d[0]
    y = cy + r * d[1]
    z = cz + r * d[2]
    return x, y, z
'''