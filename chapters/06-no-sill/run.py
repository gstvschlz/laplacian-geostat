"""Letting κ go to zero: the sill leaves and a variogram stays."""

import matplotlib.pyplot as plt
import numpy as np

import lapgeo as lg
from lapgeo import style

style.use()
kappas = [0.2, 0.05, 0.0125, 0.0]
shades = style.SEQUENTIAL(np.linspace(0.35, 1.0, len(kappas)))

fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.5))
ax = axes[0]
h = np.arange(1, 201)
for kappa, shade in zip(kappas, shades):
    c = lg.covariance((8192,), kappa, 1)
    ax.plot(h, c[0] - c[h], color=shade, label=f"κ = {kappa:g}")
ax.set_title("1-D, α = 1: the limit is the linear variogram h/2")
ax.set_xlabel("lag h, cells")
ax.set_ylabel("variogram")
ax.legend()

ax = axes[1]
h = np.unique(np.round(np.geomspace(1, 400, 40)).astype(int))
for kappa, shade in zip(kappas, shades):
    c = lg.covariance((2048, 2048), kappa, 1)
    ax.semilogx(h, c[0, 0] - c[0, h], color=shade, label=f"κ = {kappa:g}")
ax.set_title("2-D, α = 1: the limit is logarithmic (de Wijsian)")
ax.set_xlabel("lag h, cells (log scale)")
ax.legend()
fig.tight_layout()
style.save(fig, __file__, "no-sill.png")
