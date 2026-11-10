"""Unconditional draw by one sparse solve, kriging mean, and conditional draws."""

import matplotlib.pyplot as plt
import numpy as np

import lapgeo as lg
from lapgeo import style

style.use()
shape = (120, 120)
kappa = 0.12
tau = 1 / lg.matern_variance(kappa, 2, 2)
rng = np.random.default_rng(11)

truth = lg.simulate_spde(shape, kappa, tau, rng)
data = rng.choice(truth.size, 40, replace=False)
values = truth.ravel()[data]
q = lg.precision(shape, kappa, 2, tau)

mean = lg.conditional_mean(q, data, values).reshape(shape)
draws = [
    lg.conditional_simulation(q, data, values, lg.simulate_spde(shape, kappa, tau, rng)).reshape(shape)
    for _ in range(2)
]

limit = float(np.abs(truth).max())
iy, ix = np.unravel_index(data, shape)
fig, axes = plt.subplots(1, 4, figsize=(11, 3.1))
titles = ["Reference draw", "Kriging mean", "Conditional draw 1", "Conditional draw 2"]
for ax, field, title in zip(axes, [truth, mean, *draws], titles):
    style.image(ax, field, title, cmap=style.DIVERGING, vmin=-limit, vmax=limit)
    ax.scatter(
        ix,
        iy,
        c=values,
        cmap=style.DIVERGING,
        vmin=-limit,
        vmax=limit,
        s=16,
        edgecolors=style.INK,
        linewidths=0.6,
    )
fig.tight_layout()
style.save(fig, __file__, "simulation.png")
