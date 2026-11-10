"""Changing the links: anisotropy, a varying range, and a barrier."""

import matplotlib.pyplot as plt
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu

import lapgeo as lg
from lapgeo import style

style.use()
shape = (120, 120)
n = shape[0] * shape[1]
rng = np.random.default_rng(2)

aniso = lg.simulate_spde(shape, 0.25, rng=rng, weights=(1.0, 9.0))

kappa = np.broadcast_to(np.geomspace(0.6, 0.06, shape[1]), shape).ravel()
varying = splu((sp.diags_array(kappa**2) + lg.laplacian(shape, periodic=False)).tocsc()).solve(
    kappa * rng.standard_normal(n)
)

edges = lg.grid_edges(shape)
iy, ix = np.unravel_index(edges, shape)
crosses = (ix.min(axis=1) == 59) & (ix.max(axis=1) == 60) & (iy[:, 0] >= 40)
q = lg.precision(kappa=0.06, alpha=2, lap=lg.laplacian_from_edges(n, edges[~crosses]))
source = np.ravel_multi_index((55, 52), shape)
unit = np.zeros(n)
unit[source] = 1.0
column = splu(q.tocsc()).solve(unit)

fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6))
style.image(axes[0], aniso, "Links 9× stronger along x", diverging=True)
style.image(axes[1], varying.reshape(shape), "κ falling from left to right", diverging=True)
style.image(axes[2], column.reshape(shape) / column.max(), "Covariance with one cell, across a fault")
axes[2].plot([59.5, 59.5], [39.5, 119.5], color=style.INK, lw=1.5)
axes[2].plot(52, 55, "o", color=style.ORANGE, mec=style.SURFACE, mew=1.2, ms=7)
fig.tight_layout()
style.save(fig, __file__, "anisotropy-and-barriers.png")
