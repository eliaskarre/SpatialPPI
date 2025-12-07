def spring_layout_3d_constrained(
    G,
    constraint,
    iterations: int = 350,
    seed: int = 1,
    k: Optional[float] = None,
    threshold: float = 1e-4,
    recenter_each_iter: bool = True,
    center=(0.0, 0.0, 0.0),
):
    """Fruchterman-Reingold 3D + constraint forces"""
    nodes = list(G)
    n = len(nodes)
    if n == 0:
        return {}

    A = nx.to_numpy_array(G, nodelist=nodes, weight=None, dtype=float)
    if k is None:
        k = math.sqrt(1.0 / n)

    rng = np.random.default_rng(seed)
    pos = constraint.sample(n, rng)

    # 3D-safe temperature
    t = np.ptp(pos, axis=0).max() * 0.1
    t = 0.1 if t == 0 else t
    dt = t / (iterations + 1)

    for _ in range(iterations):
        delta = pos[:, None, :] - pos[None, :, :]
        dist = np.linalg.norm(delta, axis=-1)
        np.clip(dist, 0.01, None, out=dist)

        disp = np.einsum("ijk,ij->ik", delta, (k * k / dist**2 - A * dist / k))
        disp += constraint.forces(pos)

        length = np.linalg.norm(disp, axis=1)
        length = np.clip(length, 0.01, None)
        step = disp * (t / length)[:, None]
        pos += step

        if recenter_each_iter:
            pos -= pos.mean(axis=0)

        t -= dt
        if (np.linalg.norm(step) / n) < threshold:
            break

    pos += np.array(center, dtype=float)
    return {nodes[i]: pos[i] for i in range(n)}