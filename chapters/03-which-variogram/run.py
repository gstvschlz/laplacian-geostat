"""The variogram a Laplacian precision implies: Matérn with ν = α − d/2."""

import matplotlib.pyplot as plt
import numpy as np

import lapgeo as lg
from lapgeo import style

style.use()
kappa = 0.15
cases = [
    ("1-D lattice, α = 1", (4096,), 1, style.BLUE, "o"),
    ("3-D lattice, α = 2", (128, 128, 128), 2, style.AQUA, "s"),
    ("2-D lattice, α = 2", (512, 512), 2, style.ORANGE, "o"),
]

fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), width_ratios=[1.3, 1, 1])
ax = axes[0]
h = np.arange(0, 41)
fine = np.linspace(0, 40, 300)
ax.plot(fine, 1 - lg.matern(fine, kappa, 0.5), color=style.BLUE, label="exponential (ν = ½)")
ax.plot(fine, 1 - lg.matern(fine, kappa, 1.0), color=style.ORANGE, label="Matérn ν = 1")
for label, shape, alpha, colour, marker in cases:
    d = len(shape)
    sill = lg.matern_variance(kappa, alpha, d)
    lattice = lg.covariance(shape, kappa, alpha)[(0,) * (d - 1)][h]
    pick = slice(1, None, 4) if marker == "s" else slice(3, None, 4)
    ax.plot(
        h[pick], 1 - lattice[pick] / sill, marker, color=colour, mec=style.SURFACE, mew=1, ms=6, label=label
    )
ax.set_title("Lattice (markers) against Matérn (lines)")
ax.set_xlabel("lag h, cells")
ax.set_ylabel("variogram / sill")
ax.set_ylim(0)
ax.legend(loc="lower right")

shape = (160, 160)
for ax, alpha in zip(axes[1:], (1, 2)):
    z = lg.simulate_spectral(shape, kappa, alpha, rng=3)
    style.image(ax, z, f"2-D draw, α = {alpha} (ν = {alpha - 1})", diverging=True)
fig.tight_layout()
style.save(fig, __file__, "which-variogram.png")
