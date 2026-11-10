"""Covariance is dense, precision is the stencil, kriging weights are the stencil."""

import matplotlib.pyplot as plt
import numpy as np

import lapgeo as lg
from lapgeo import style

style.use()
kappa = 0.3

n = 400
c1 = lg.covariance((n,), kappa, 1)
rho = 1 + kappa**2 / 2 - np.sqrt(kappa**2 + kappa**4 / 4)
h = np.arange(0, 21)

shape = (21, 21)
centre = np.ravel_multi_index((10, 10), shape)
q = lg.precision(shape, kappa, 1).toarray()
cov = np.linalg.inv(q)
others = np.delete(np.arange(q.shape[0]), centre)
_, w = lg.simple_kriging(cov, others, np.zeros(others.size))
weights = np.zeros(q.shape[0])
weights[others] = w[centre]

fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), width_ratios=[1.25, 1, 1])
ax = axes[0]
ax.plot(h, c1[0] * rho**h, color=style.BLUE, label="exponential, C(0) ρ^h")
ax.plot(h, c1[h], "o", color=style.ORANGE, mec=style.SURFACE, mew=1, label="lattice")
ax.set_title("1-D: inverse of κ²I + L")
ax.set_xlabel("lag h, cells")
ax.set_ylabel("covariance")
ax.set_ylim(0)
ax.legend()

style.image(axes[1], cov[centre].reshape(shape), "2-D: covariance with a cell")
im = style.image(axes[2], weights.reshape(shape), "2-D: kriging weights for it", diverging=True)
axes[2].text(10, 12.2, f"{weights.max():.3f} each", ha="center", color=style.INK)
fig.tight_layout()
style.save(fig, __file__, "screen-effect.png")
