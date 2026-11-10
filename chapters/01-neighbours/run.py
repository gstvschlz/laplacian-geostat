"""A surface, its Laplacian, and the Laplacian's eigenvalues."""

import matplotlib.pyplot as plt
import numpy as np

import lapgeo as lg
from lapgeo import style

style.use()
shape = (96, 96)
y, x = np.indices(shape)


def bump(cy, cx, width):
    return np.exp(-((y - cy) ** 2 + (x - cx) ** 2) / (2 * width**2))


z = 0.012 * x + bump(62, 30, 9) - 0.8 * bump(30, 66, 12)
lz = (lg.laplacian(shape, periodic=False) @ z.ravel()).reshape(shape)

fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), width_ratios=[1, 1, 1.15])
style.image(axes[0], z[1:-1, 1:-1], "z: a ramp, a high and a low", diverging=True)
style.image(axes[1], lz[1:-1, 1:-1], "Lz = 4z − sum of neighbours", diverging=True)
axes[1].text(29, 61, "+", ha="center", va="center", color=style.SURFACE, fontsize=13)
axes[1].text(65, 29, "−", ha="center", va="center", color=style.SURFACE, fontsize=13)
axes[1].text(47, 86, "zero on the ramp", ha="center", color=style.MUTED)

omega = np.linspace(0, np.pi, 200)
ax = axes[2]
ax.plot(omega, omega**2, color=style.MUTED, lw=1.2, ls="--")
ax.plot(omega, 2 - 2 * np.cos(omega), color=style.BLUE)
ax.text(2.25, 6.3, "ω²  (continuous −Δ)", color=style.MUTED, ha="right")
ax.text(3.1, 1.7, "2 − 2 cos ω  (lattice L)", color=style.INK, ha="right")
ax.set_title("What L does to a wave")
ax.set_xlabel("ω, radians per cell (π = a wave of two cells)")
ax.set_ylabel("factor multiplying the wave")
ax.set_ylim(0, 10)
fig.tight_layout()
style.save(fig, __file__, "neighbours.png")
