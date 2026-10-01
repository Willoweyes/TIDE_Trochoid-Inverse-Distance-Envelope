# -*- coding: utf-8 -*-
"""Shared publication style for every figure.

EDIT THIS FILE to restyle all figures at once: colours, fonts, line widths,
marker set, figure sizes and DPI are all defined here and nowhere else.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- palette (colour-blind safe, prints legibly in greyscale) -------------
C_EXACT    = "#1f4e79"   # TIDE / exact reference
C_FORWARD  = "#c1440e"   # forward Z-map
C_THIRD    = "#2e7d32"
C_FOURTH   = "#7b5aa6"
C_GREY     = "#5a5a5a"
C_LIGHT    = "#b0b0b0"
CYCLE      = [C_EXACT, C_FORWARD, C_THIRD, C_FOURTH, C_GREY]

# ---- sizes (inches): Springer LNCS/LNME text width is 122 mm = 4.80 in ------
W_FULL = 4.80
W_ONE  = 3.4      # (legacy two-column layouts)
W_TWO  = 7.0
DPI    = 600

MARKERS = ["o", "s", "^", "D", "v", "P", "X"]

RC = {
    "figure.dpi":            120,
    "savefig.dpi":           DPI,
    "savefig.bbox":          "tight",
    "savefig.pad_inches":    0.02,
    "font.family":           "serif",
    "font.serif":            ["TeX Gyre Termes", "STIXGeneral", "Times New Roman", "DejaVu Serif"],
    "font.size":             8.0,
    "axes.titlesize":        8.0,
    "axes.labelsize":        8.0,
    "xtick.labelsize":       7.0,
    "ytick.labelsize":       7.0,
    "legend.fontsize":       6.8,
    "svg.fonttype":          "path",
    "pdf.fonttype":          42,
    "legend.frameon":        True,
    "legend.framealpha":     0.92,
    "legend.edgecolor":      "0.75",
    "legend.borderpad":      0.35,
    "legend.labelspacing":   0.3,
    "legend.handlelength":   1.8,
    "axes.linewidth":        0.7,
    "axes.grid":             True,
    "axes.axisbelow":        True,
    "grid.linewidth":        0.4,
    "grid.alpha":            0.35,
    "grid.color":            "0.6",
    "lines.linewidth":       1.1,
    "lines.markersize":      3.4,
    "lines.markeredgewidth": 0.6,
    "xtick.direction":       "in",
    "ytick.direction":       "in",
    "xtick.major.width":     0.7,
    "ytick.major.width":     0.7,
    "xtick.minor.width":     0.5,
    "ytick.minor.width":     0.5,
    "mathtext.fontset":      "stix",
}


def apply():
    plt.rcParams.update(RC)
    plt.rcParams["axes.prop_cycle"] = plt.cycler(color=CYCLE)
    return plt


def save(fig, root, name):
    """Write PDF, SVG (text as paths, for Word) and a 600-dpi PNG."""
    for ext in ("pdf", "svg", "png"):
        fig.savefig(root / f"paper/figures/{name}.{ext}")
