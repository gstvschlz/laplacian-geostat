"""Estimation as a differential equation between the data, in one dimension."""

import matplotlib.pyplot as plt
import numpy as np

import lapgeo as lg
from lapgeo import style

style.use()
n = 240
data = np.array([20, 55, 75, 130, 150, 205])
values = np.array([0.2, 1.0, 1.6, -0.4, -0.9, 0.5])
x = np.arange(n)

curves = [
    ("α = 1, κ = 0: straight lines", lg.precision((n,), 0.0, 1), style.BLUE),
    ("α = 1, κ = 0.08: sags to the mean", lg.precision((n,), 0.08, 1), style.AQUA),
    ("α = 2, κ = 0: cubic spline", lg.precision((n,), 0.0, 2), style.ORANGE),
]

fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.6), width_ratios=[1.6, 1])
ax = axes[0]
ax.axhline(0, color=style.MUTED, lw=0.8)
for label, q, colour in curves:
    ax.plot(x, lg.conditional_mean(q, data, values), color=colour, label=label)
ax.plot(data, values, "o", color=style.INK, mec=style.SURFACE, mew=1.2, ms=7, zorder=5, label="data")
ax.set_title("Kriging mean on a ring of 240 cells")
ax.set_xlabel("cell")
ax.legend(ncols=2, loc="lower left", fontsize=8)
ax.set_ylim(-2.6, 2.1)

shape = (64, 64)
kappa = 0.1
rng = np.random.default_rng(5)
d2 = rng.choice(shape[0] * shape[1], 30, replace=False)
v2 = rng.standard_normal(30)
exact = lg.conditional_mean(lg.precision(shape, kappa, 1), d2, v2)
sweeps = np.arange(0, 401, 20)
z = np.zeros(shape)
errors = []
for _ in sweeps:
    errors.append(np.abs(z.ravel() - exact).max() if errors else np.abs(exact).max())
    z = lg.relax(shape, kappa, d2, v2, 20, start=z)
errors[0] = np.abs(lg.relax(shape, kappa, d2, v2, 0).ravel() - exact).max()
ax = axes[1]
ax.semilogy(sweeps, errors, color=style.BLUE)
ax.set_title("Averaging neighbours converges to kriging")
ax.set_xlabel("sweeps of the cellular rule (2-D, α = 1)")
ax.set_ylabel("largest error against the solve")
fig.tight_layout()
style.save(fig, __file__, "estimation.png")
