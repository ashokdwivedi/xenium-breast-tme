"""Figure style and reusable plot helpers.

A single :func:`set_style` call fixes the typography ladder, tick direction and
output resolution for every figure in the project, so panels produced by
different scripts composite cleanly.
"""

from __future__ import annotations

import matplotlib as mpl
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Palette for the breast tumour microenvironment. Bound once here and reused
# everywhere: colour is the cross-reference between UMAP, tissue map, niche
# profile and interface plot, so it must not be re-derived per figure.
CELLTYPE_COLORS: dict[str, str] = {
    "Cancer Epithelial": "#B2182B",
    "Epithelial (tumour)": "#B2182B",
    "Normal Epithelial": "#E8A0A8",
    "T-cells": "#1F78B4",
    "T cell": "#1F78B4",
    "B-cells": "#6BAED6",
    "B cell": "#6BAED6",
    "Plasmablasts": "#9ECAE1",
    "Plasma cell": "#9ECAE1",
    "Myeloid": "#FF7F00",
    "Mast cell": "#FDBF6F",
    "Endothelial": "#33A02C",
    "CAFs": "#6A3D9A",
    "Fibroblast (CAF)": "#6A3D9A",
    "PVL": "#C8A2C8",
    "Perivascular (PVL)": "#C8A2C8",
    "Unassigned": "#BDBDBD",
}

SAMPLE_COLORS = {"Rep1": "#2B6CB0", "Rep2": "#DD6B20"}
ALARM = "#D62728"   # reserved for flags/failures; never a data-series colour


def panel_letter(ax, letter: str, dx: float = -22, dy: float = 6) -> None:
    """Bold panel letter placed outside the axes, so it never collides with a title.

    Offsets are in points rather than axes fractions: an axes-fraction offset
    scales with panel width, which is exactly wrong when panels differ in size.
    """
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction",
                xytext=(dx, dy), textcoords="offset points",
                fontweight="bold", fontsize=mpl.rcParams["font.size"] + 2,
                va="bottom", ha="left", annotation_clip=False)


def set_style(sizes: tuple[int, int, int] = (9, 8, 7), dpi: int = 300) -> None:
    """Apply the project figure style: three-size ladder, outward ticks, 300 dpi."""
    base, small, tick = sizes
    mpl.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": dpi,
        "savefig.bbox": "tight",
        "font.size": base,
        "axes.titlesize": base,
        "axes.labelsize": base,
        "axes.titlelocation": "left",
        "axes.titleweight": "regular",
        "legend.fontsize": small,
        "legend.frameon": False,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def celltype_palette(labels) -> dict[str, str]:
    """Colour map for a set of labels, falling back to tab20 for unknown names."""
    labels = list(pd.unique(pd.Series(labels).astype(str)))
    fallback = plt.get_cmap("tab20").colors
    out, spare = {}, 0
    for lab in labels:
        if lab in CELLTYPE_COLORS:
            out[lab] = CELLTYPE_COLORS[lab]
        else:
            out[lab] = mpl.colors.to_hex(fallback[spare % len(fallback)])
            spare += 1
    return out


def spatial_scatter(ax, xy: np.ndarray, values, *, categorical: bool,
                    palette: dict | None = None, cmap: str = "viridis",
                    s: float = 0.6, vmin=None, vmax=None, rasterized: bool = True):
    """Tissue map. Equal aspect, no ticks, a scale bar instead of axis numbers."""
    if categorical:
        values = pd.Series(values).astype(str).to_numpy()
        palette = palette or celltype_palette(values)
        for lab, colour in palette.items():
            m = values == lab
            if m.any():
                ax.scatter(xy[m, 0], xy[m, 1], s=s, c=colour, lw=0,
                           rasterized=rasterized, label=lab)
        handle = None
    else:
        handle = ax.scatter(xy[:, 0], xy[:, 1], s=s, c=values, cmap=cmap, lw=0,
                            vmin=vmin, vmax=vmax, rasterized=rasterized)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.invert_yaxis()  # image convention: y increases downward on the slide
    return handle


def add_scalebar(ax, xy: np.ndarray, length_um: float = 1000, label: str | None = None) -> None:
    """Horizontal scale bar in the lower-left of a tissue map."""
    x0, x1 = xy[:, 0].min(), xy[:, 0].max()
    y0, y1 = xy[:, 1].min(), xy[:, 1].max()
    pad_x, pad_y = 0.04 * (x1 - x0), 0.06 * (y1 - y0)
    xs = x0 + pad_x
    ys = y1 - pad_y
    # White halo: tissue maps have both dark and light backgrounds, and the bar
    # has to stay legible over either without knowing which it lands on.
    halo = [pe.withStroke(linewidth=2.2, foreground="white")]
    ax.plot([xs, xs + length_um], [ys, ys], color="black", lw=1.6,
            solid_capstyle="butt", path_effects=halo)
    ax.text(xs + length_um / 2, ys - 0.012 * (y1 - y0),
            label or f"{length_um:g} \u00b5m", ha="center", va="bottom",
            fontsize=plt.rcParams["legend.fontsize"], path_effects=halo)


def bbox_report(fig) -> list[tuple[str, str]]:
    """Geometric overlap check: returns colliding visible text pairs."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible()]
    return [(a.get_text(), b.get_text())
            for i, (a, ba) in enumerate(texts) for b, bb in texts[i + 1:] if ba.overlaps(bb)]
