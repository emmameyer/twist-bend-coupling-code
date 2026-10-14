"""
fig_style.py

script to ensure same format across all figures
import in every plotter to keep marker sizes, linewidths, fontsizes, tick geometry, etc. consistent


to use, put this in imports:
    from fig_style import RCPARAMS, MS, LW, apply_style, FIG_SIZE_MM

    plt.rcParams.update(RCPARAMS)
    # or you could do:
    apply_style()
"""

import matplotlib as mpl
import matplotlib.pyplot as plt

# figure size

FIG_SIZE_MM = (53, 37)          # (width, height) in millimetres, specific measurement needed for some paper figures


def fig_size_inches(width_mm=FIG_SIZE_MM[0], height_mm=FIG_SIZE_MM[1]):
    """Return (width, height) in inches."""
    return width_mm / 25.4, height_mm / 25.4


# marker sizes

class MS:
    """Marker size constants."""
    small  = 1.5
    medium = 2.5    
    large  = 3.5
    cross  = 3.5


# linewidths

class LW:
    """Line width constants."""
    axes   = 0.5   # spine / frame
    tick   = 0.5   # tick marks
    data   = 1.0   # data lines (curves, parabola segments)
    split  = 1.00   # vertical split bars that we used in the sawtooth plots
    grid   = 0.75   # background grid / reference lines
    edge   = 0.75   # marker edge (used mostly for the hollow points)
    fit    = 1.00   # polynomial fit line


# ticks

class TICK:
    major_length = 3.5
    minor_length = 2.0
    pad_major    = 1.0
    pad_minor    = 1.0


# rcParams dict

RCPARAMS = {
    # fonts
    "svg.fonttype":    "path",
    "font.family":          "sans-serif",
    "font.sans-serif":      "Arial",
    "font.size":             6,
    "axes.labelsize":        6,
    "axes.titlesize":        6,
    "xtick.labelsize":       6,
    "ytick.labelsize":       6,
    # axes
    "axes.linewidth":        LW.axes,
    # ticks
    "xtick.major.pad":       TICK.pad_major,
    "ytick.major.pad":       TICK.pad_major,
    "xtick.minor.pad":       TICK.pad_minor,
    "ytick.minor.pad":       TICK.pad_minor,
    "xtick.major.size":      TICK.major_length,
    "ytick.major.size":      TICK.major_length,
    "xtick.minor.size":      TICK.minor_length,
    "ytick.minor.size":      TICK.minor_length,
    "xtick.major.width":     LW.tick,
    "ytick.major.width":     LW.tick,
    "xtick.minor.width":     LW.tick,
    "ytick.minor.width":     LW.tick,
    "xtick.direction":       "out",
    "ytick.direction":       "out",
    # label pads
    "axes.labelpad":         2.0,
}


def apply_style():
    """Apply RCPARAMS to the current matplotlib session."""
    plt.rcParams.update(RCPARAMS)
