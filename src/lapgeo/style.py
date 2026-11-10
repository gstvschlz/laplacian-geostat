"""Figure style shared by the chapters."""

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

SEQUENTIAL = LinearSegmentedColormap.from_list("lap_seq", ["#f3f6fb", "#a9c8ef", BLUE, "#123a6b"])
DIVERGING = LinearSegmentedColormap.from_list("lap_div", ["#123a6b", BLUE, "#ececea", ORANGE, "#7a2c0e"])


def use():
    mpl.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "savefig.dpi": 150,
            "savefig.bbox": "tight",
            "font.size": 9,
            "axes.titlesize": 9.5,
            "axes.titleweight": "bold",
            "axes.titlelocation": "left",
            "axes.edgecolor": MUTED,
            "axes.labelcolor": MUTED,
            "axes.linewidth": 0.6,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.color": GRID,
            "grid.linewidth": 0.6,
            "axes.axisbelow": True,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "text.color": INK,
            "lines.linewidth": 2.0,
            "lines.markersize": 5,
            "legend.frameon": False,
            "image.cmap": "lap_seq",
            "image.interpolation": "nearest",
            "axes.prop_cycle": mpl.cycler(color=[BLUE, ORANGE, AQUA]),
        }
    )
    for cmap in (SEQUENTIAL, DIVERGING):
        if cmap.name not in mpl.colormaps:
            mpl.colormaps.register(cmap)


def image(ax, values, title, diverging=False, **kwargs):
    """Show a grid without axes; a diverging map is centred on zero."""
    if diverging:
        limit = float(abs(values).max())
        kwargs = {"cmap": DIVERGING, "vmin": -limit, "vmax": limit} | kwargs
    im = ax.imshow(values, origin="lower", **kwargs)
    ax.set_title(title)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for side in ax.spines.values():
        side.set_visible(False)
    return im


def save(fig, script, name):
    path = Path(script).resolve().parent / name
    fig.savefig(path)
    plt.close(fig)
    print(path.relative_to(Path(script).resolve().parents[2]))
