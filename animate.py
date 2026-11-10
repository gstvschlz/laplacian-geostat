"""The Laplacian in geostatistical estimation and simulation, as one animation.

Writes laplacian.mp4 beside this file.
"""

from itertools import pairwise
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FFMpegWriter

import lapgeo as lg
from lapgeo import style

FPS = 20
HERE = Path(__file__).resolve().parent
style.use()
plt.rcParams.update(
    {"font.size": 13, "axes.titlesize": 13, "axes.titleweight": "normal", "savefig.bbox": None}
)


def ease(t):
    t = np.clip(t, 0.0, 1.0)
    return t * t * (3 - 2 * t)


def schedule(frames, last):
    return np.unique(np.r_[0, np.round(np.geomspace(1, last, frames - 1)).astype(int)])


def captions(fig, title, sub, step):
    fig.text(0.05, 0.93, title, fontsize=20, fontweight="bold", va="center")
    fig.text(0.05, 0.865, sub, fontsize=14.5, color=style.MUTED, va="center")
    fig.text(0.985, 0.025, f"{step} / 6", fontsize=12, color=style.MUTED, va="center", ha="right")


def footer(fig, text):
    fig.text(0.5, 0.04, text, fontsize=14, ha="center", va="center")


def second_difference(z):
    out = np.zeros_like(z)
    out[1:-1] = 2 * z[1:-1] - z[:-2] - z[2:]
    return out


# ---------------------------------------------------------------- 1-D scenes

N = 49
X = np.arange(N)
PROFILE = 0.02 * X + 1.3 * np.exp(-(((X - 16) / 2.6) ** 2)) - 1.0 * np.exp(-(((X - 34) / 3.0) ** 2))
DATA = np.array([4, 14, 22, 33, 44])
VALUES = np.array([0.4, 1.3, 0.6, -0.9, 0.3])
FREE = np.ones(N, dtype=bool)
FREE[DATA] = False


def profile_axes(fig, lower):
    if lower:
        top = fig.add_axes([0.07, 0.36, 0.86, 0.43])
        bottom = fig.add_axes([0.07, 0.09, 0.86, 0.19], sharex=top)
        bottom.set_ylim(-0.42, 0.42)
        bottom.set_ylabel("Laplacian")
        bottom.axhline(0, color=style.MUTED, lw=0.8)
        plt.setp(top.get_xticklabels(), visible=False)
        bottom.set_xlabel("cell")
    else:
        top = fig.add_axes([0.07, 0.11, 0.86, 0.68])
        top.set_xlabel("cell")
        bottom = None
    top.set_xlim(-1, N)
    top.set_ylim(-1.35, 2.25)
    top.set_ylabel("value")
    return top, bottom


def bars(ax, x, values):
    colours = np.where(values >= 0, style.ORANGE, style.BLUE)
    ax.bar(x, values, width=0.7, color=colours)


def scene_operator(fig, t):
    captions(
        fig,
        "The Laplacian of a cell: its value minus the average of its neighbours",
        "Positive on a high, negative in a low, zero on a flat or a straight slope.",
        1,
    )
    top, bottom = profile_axes(fig, True)
    i = int(np.clip(1 + t / 0.92 * (N - 2), 1, N - 2))
    lap = second_difference(PROFILE)
    top.plot(X, PROFILE, color=style.GRID, lw=1.5, zorder=1)
    top.plot(X, PROFILE, "o", color=style.MUTED, ms=5, zorder=2)
    middle = 0.5 * (PROFILE[i - 1] + PROFILE[i + 1])
    top.plot([i - 1, i + 1], [PROFILE[i - 1], PROFILE[i + 1]], color=style.BLUE, lw=2, zorder=3)
    top.plot([i - 1, i + 1], [PROFILE[i - 1], PROFILE[i + 1]], "o", color=style.BLUE, ms=10, zorder=4)
    top.plot([i, i], [middle, PROFILE[i]], color=style.ORANGE, lw=4, zorder=3, solid_capstyle="butt")
    top.plot(i, middle, "o", mfc=style.SURFACE, mec=style.BLUE, mew=2, ms=9, zorder=5)
    top.plot(i, PROFILE[i], "o", color=style.ORANGE, mec=style.SURFACE, mew=1.5, ms=12, zorder=6)
    top.text(
        0.99,
        0.95,
        f"value {PROFILE[i]:+.2f}    average of neighbours {middle:+.2f}    difference {PROFILE[i] - middle:+.3f}",
        transform=top.transAxes,
        ha="right",
        va="top",
        fontsize=13,
    )
    bars(bottom, X[1 : i + 1], lap[1 : i + 1])


RELAX_SWEEPS = schedule(110, 150)


def relax_line():
    z = np.zeros(N)
    z[DATA] = VALUES
    states = {0: z.copy()}
    for sweep in range(1, RELAX_SWEEPS[-1] + 1):
        for parity in (0, 1):
            new = z.copy()
            new[1:-1] = 0.5 * (z[:-2] + z[2:])
            new[0], new[-1] = z[1], z[-2]
            pick = FREE & (X % 2 == parity)
            z[pick] = new[pick]
        if sweep in RELAX_SWEEPS:
            states[sweep] = z.copy()
    return [states[s] for s in RELAX_SWEEPS]


LINE_STATES = relax_line()


def scene_estimate(fig, t):
    captions(
        fig,
        "Estimation: hold the data, make the Laplacian zero everywhere else",
        "Rule: every cell without a datum becomes the average of its neighbours. Repeat.",
        2,
    )
    top, bottom = profile_axes(fig, True)
    k = min(int(t / 0.85 * len(LINE_STATES)), len(LINE_STATES) - 1)
    z = LINE_STATES[k]
    top.axhline(0, color=style.MUTED, lw=0.8)
    top.plot(X, z, color=style.BLUE, lw=2, zorder=2)
    top.plot(X[FREE], z[FREE], "o", color=style.BLUE, ms=5, zorder=3)
    top.plot(DATA, VALUES, "o", color=style.INK, mec=style.SURFACE, mew=1.5, ms=11, zorder=5)
    done = k == len(LINE_STATES) - 1
    note = "straight lines between the data: kriging with a linear variogram" if done else ""
    top.text(0.99, 0.95, f"sweep {RELAX_SWEEPS[k]}", transform=top.transAxes, ha="right", va="top")
    top.text(0.01, 0.95, note, transform=top.transAxes, ha="left", va="top", color=style.INK)
    lap = second_difference(z)
    lap[~FREE] = 0
    bars(bottom, X[1:-1], np.clip(lap[1:-1], -0.42, 0.42))


LINE_LAPLACIAN = lg.laplacian((N,), periodic=False)


def scene_kappa(fig, t):
    captions(
        fig,
        "Add κ²: each cell becomes a little less than the average",
        "The estimate sags to the mean away from the data. κ sets the range: simple kriging, exponential variogram.",
        3,
    )
    top, _ = profile_axes(fig, False)
    kappa = 0.4 * ease(t / 0.8)
    z = lg.conditional_mean(lg.precision(kappa=kappa, alpha=1, lap=LINE_LAPLACIAN), DATA, VALUES)
    top.axhline(0, color=style.MUTED, lw=0.8)
    top.text(N - 0.5, 0.06, "mean", color=style.MUTED, ha="right", fontsize=11)
    top.plot(X, LINE_STATES[-1], color=style.GRID, lw=2, zorder=1)
    top.plot(X, z, color=style.BLUE, lw=2.5, zorder=2)
    top.plot(DATA, VALUES, "o", color=style.INK, mec=style.SURFACE, mew=1.5, ms=11, zorder=5)
    reach = "no range, no sill" if kappa < 0.02 else f"practical range ≈ 3/κ = {3 / kappa:.0f} cells"
    top.text(0.99, 0.95, f"κ = {kappa:.2f}    {reach}", transform=top.transAxes, ha="right", va="top")
    top.text(
        0.01,
        0.95,
        f"rule: cell = sum of 2 neighbours / (2 + κ²) = {2 / (2 + kappa**2):.3f} × their average",
        transform=top.transAxes,
        ha="left",
        va="top",
    )


# ---------------------------------------------------------------- 2-D scenes

SHAPE = (96, 96)
RNG = np.random.default_rng(12)
KAPPA = 0.12
TAU = lg.matern_variance(KAPPA, 2, 2)
Q2 = lg.precision(SHAPE, KAPPA, 2, TAU)
REFERENCE = lg.simulate_spectral(SHAPE, KAPPA, 2, TAU, RNG)
DATA2 = RNG.choice(REFERENCE.size, 30, replace=False)
VALUES2 = REFERENCE.ravel()[DATA2]
DY, DX = np.unravel_index(DATA2, SHAPE)
LIMIT = 2.6

MAP_SWEEPS = schedule(130, 2500)


def relax_map(source, data, values, sweeps, kappa):
    z = np.zeros(SHAPE)
    states = [lg.relax(SHAPE, kappa, data, values, 0, start=z)]
    for a, b in pairwise(sweeps):
        states.append(lg.relax(SHAPE, kappa, data, values, b - a, start=states[-1], source=source))
    return states


KRIGING_KAPPA = 0.05
MAP_STATES = relax_map(None, DATA2, VALUES2, MAP_SWEEPS, KRIGING_KAPPA)

NOISE = RNG.standard_normal(SHAPE)
NOISE_KAPPA = 0.15
NOISE_SWEEPS = schedule(140, 1500)
NOISE_STATES = relax_map(NOISE, [], [], NOISE_SWEEPS, NOISE_KAPPA)

MEAN2 = lg.conditional_mean(Q2, DATA2, VALUES2).reshape(SHAPE)
RESIDUALS = []
for _ in range(4):
    draw = lg.simulate_spectral(SHAPE, KAPPA, 2, TAU, RNG)
    RESIDUALS.append(draw - lg.conditional_mean(Q2, DATA2, draw.ravel()[DATA2]).reshape(SHAPE))


def panel(fig, box, values, title, data=True):
    ax = fig.add_axes(box)
    style.image(ax, values, title, cmap=style.DIVERGING, vmin=-LIMIT, vmax=LIMIT)
    if data:
        ax.scatter(
            DX,
            DY,
            c=VALUES2,
            cmap=style.DIVERGING,
            vmin=-LIMIT,
            vmax=LIMIT,
            s=70,
            edgecolors=style.INK,
            linewidths=1.2,
        )
    return ax


LEFT, RIGHT = [0.07, 0.09, 0.40, 0.71], [0.53, 0.09, 0.40, 0.71]


def stencil(fig):
    ax = fig.add_axes(LEFT)
    ax.set_xlim(-2.6, 2.6)
    ax.set_ylim(-2.6, 2.6)
    ax.set_aspect("equal")
    ax.axis("off")
    for dy, dx, label in [(0, 0, "4"), (0, 1, "−1"), (0, -1, "−1"), (1, 0, "−1"), (-1, 0, "−1")]:
        centre = dy == 0 and dx == 0
        colour = style.ORANGE if centre else style.BLUE
        ax.add_patch(plt.Rectangle((dx - 0.46, dy - 0.46), 0.92, 0.92, color=colour))
        ax.text(dx, dy, label, ha="center", va="center", color=style.SURFACE, fontsize=22, fontweight="bold")
    ax.text(0, -2.15, "L: four times the cell minus its four neighbours", ha="center", fontsize=14)
    ax.text(0, 2.05, "the same operator on a map", ha="center", fontsize=14, color=style.MUTED)


def scene_map(fig, t):
    captions(
        fig,
        "On a map, the same rule gives the kriging map",
        "Every cell without a datum becomes the sum of its four neighbours / (4 + κ²). Repeat.",
        4,
    )
    stencil(fig)
    k = min(int(t / 0.85 * len(MAP_STATES)), len(MAP_STATES) - 1)
    panel(fig, RIGHT, MAP_STATES[k], f"estimate after {MAP_SWEEPS[k]} sweeps (operator applied once, α = 1)")


def scene_simulate(fig, t):
    captions(
        fig,
        "Simulation: the same operator, fed with white noise",
        "Solve (κ² I + L) z = noise. Each cell: (sum of neighbours + its own random number) / (4 + κ²). Repeat.",
        5,
    )
    k = min(int(t / 0.85 * len(NOISE_STATES)), len(NOISE_STATES) - 1)
    z = NOISE_STATES[k]
    panel(fig, LEFT, NOISE, "white noise, one number per cell", data=False)
    panel(fig, RIGHT, z / (z.std() if k else 1.0), f"z after {NOISE_SWEEPS[k]} sweeps", data=False)
    footer(fig, f"The result has a Matérn variogram (ν = 1) with range ≈ 4/κ = {4 / NOISE_KAPPA:.0f} cells.")


def scene_conditional(fig, t):
    captions(
        fig,
        "Both at once: hold the data and feed noise",
        "Without noise it is kriging. With noise, every frame is a conditional simulation that honours the data.",
        6,
    )
    s = min(t, 0.999) * len(RESIDUALS)
    a = int(s)
    b = (a + 1) % len(RESIDUALS)
    angle = 0.5 * np.pi * ease(s - a)
    field = MEAN2 + np.cos(angle) * RESIDUALS[a] + np.sin(angle) * RESIDUALS[b]
    panel(fig, LEFT, MEAN2, "no noise: kriging (operator applied twice, α = 2)")
    panel(fig, RIGHT, field, "noise: conditional simulation (α = 2)")
    footer(fig, "precision  Q = τ (κ² I + L)^α        κ: range        α: smoothness        τ: sill")


SCENES = [
    (scene_operator, 10),
    (scene_estimate, 9),
    (scene_kappa, 8),
    (scene_map, 9),
    (scene_simulate, 9),
    (scene_conditional, 12),
]


def render(path=HERE / "laplacian.mp4", only=None):
    fig = plt.figure(figsize=(12.8, 7.2), dpi=100)
    writer = FFMpegWriter(fps=FPS, codec="libx264", extra_args=["-pix_fmt", "yuv420p", "-crf", "20"])
    with writer.saving(fig, str(path), dpi=100):
        for number, (scene, seconds) in enumerate(SCENES):
            if only is not None and number not in only:
                continue
            frames = seconds * FPS
            for frame in range(frames):
                fig.clear()
                scene(fig, frame / (frames - 1))
                writer.grab_frame(facecolor=style.SURFACE)
    plt.close(fig)
    return path


def still(scene, t, path):
    fig = plt.figure(figsize=(12.8, 7.2), dpi=100)
    scene(fig, t)
    fig.savefig(path, dpi=100)
    plt.close(fig)


if __name__ == "__main__":
    print(render().name)
