"""
plot_combo_paper.py

plots sawtooth plots for a variety of Soumya's experimental BAL-31 digest data: no additives, 5 ug/mL EtBr, 12.5 ug/mL EtBr, and HMfB
uses txt files as input (bal31 holds no additives and etbr data, the hmfb is in a separate text file)

Usage:
    python3 plot_combo_paper.py output --bal31-data bal31_data.txt --hmfb-data hmfb_data.txt --xlim 50,106
"""

import argparse
import os
import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import defaultdict
from matplotlib.lines import Line2D
from matplotlib.ticker import NullFormatter, StrMethodFormatter
from fig_style import RCPARAMS, MS, LW, TICK, apply_style, fig_size_inches
apply_style()

# constants
h_canonical  = 10.5   # canonical helical repeat (bp/turn)
rise_per_bp  = 0.34   # nm per base pair for B-DNA
default_y_lim = (9.7, 14.5)

colors = {
    "fit_line": "#4d4dff", "light_green_split": "green", "dark_green_split": "#035803", "green_split": "green", "theory_dots": "lightgray",
    "etbr_band": "blue", "etbr_light_split": "#58b4ffff", "etbr_dark_split": "#0055ff", "hmfb_band": "red", "hmfb_light_split": "#fa9595",
}

condition_label_colors = {"No Additives": "green", "EtBr 5 ug/mL": "blue", "EtBr 12.5 ug/mL": "blue", "HMfB": "red"}

custom_splits      = {"before": {63, 83, 93}, "after": {63, 83, 93}}
green_splits       = {73, 83, 92, 93, 103}
condition_splits   = {"EtBr 5 ug/mL": [54, 64, 74, 84, 94, 103], "EtBr 12.5 ug/mL": [54, 65, 75, 87, 98], "HMfB": [71, 81, 91, 102]}
condition_split_exclusions = {"EtBr 5 ug/mL": {102, 103, 104, 105}}

faint_points = {
    (57,  4,  "EtBr 5 ug/mL",    "after"),
    (63,  6,  "EtBr 5 ug/mL",    "after"),
    (92,  8,  "EtBr 5 ug/mL",    "after"),
    (55,  5,  "EtBr 12.5 ug/mL", "after"),
    (64,  6,  "EtBr 12.5 ug/mL", "after"),
    (74,  6,  "EtBr 12.5 ug/mL", "after"),
    (102, 9,  "EtBr 5 ug/mL",    "after"),
    (102, 10, "EtBr 5 ug/mL",    "after"),
    (103, 9,  "EtBr 5 ug/mL",    "after"),
    (103, 10, "EtBr 5 ug/mL",    "after"),
    (104, 9,  "EtBr 5 ug/mL",    "after"),
    (104, 10, "EtBr 5 ug/mL",    "after"),
    (105, 9,  "EtBr 5 ug/mL",    "after"),
    (105, 10, "EtBr 5 ug/mL",    "after"),
    (52,  4,  "HMfB",            "before"),
    (53,  5,  "HMfB",            "before"),
    (71,  6,  "HMfB",            "before"),
    (81,  7,  "HMfB",            "before"),
    (91,  9,  "HMfB",            "before"),
    (65,  6,  "HMfB",            "after"),
    (76,  7,  "HMfB",            "after"),
}

before_arrows = {
    "12.5": [((88, 11.3), (82, 11.3), "blue", 1.8), ((95, 11.75), (95, 10.6), "blue", 1.8)],
    "HMfB": [((80, 11.0), (83.5, 11.0), "red", 1.3), ((72, 10.5), (72, 12.0), "red", 1.1)]}

# data parsing
def parse_bal31_data(path: str) -> dict:
    """
    Parse the no adds & etbr bal-31 digest data file.
    Returns a dict with keys:
        "initial_population": list of (n_bp, lk)
        "turning_points":     list of n_bp (No Additives only)
        "conditions":         dict of condition_label -> {before, after, turning_points}
    """
    data = {"initial_population": [], "turning_points": [], "conditions": {}}
    current_section_key = None
    current_state = None

    with open(path, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.lower().startswith("# theory"):
                continue

            if line.startswith("# N_bp") and "h_expected" in line:
                current_section_key = "initial_population"
                current_state = None
                continue

            match = re.match(r"# (.*?) (before digest|after digest|turning points)", line, re.IGNORECASE)
            if match:
                label = match.group(1).strip()
                current_state = match.group(2).strip().lower()
                if label.lower() in ("no additives", "no additive"):
                    label = "No Additives"
                if label not in data["conditions"]:
                    data["conditions"][label] = {"before": [], "after": [], "turning_points": []}
                current_section_key = label
                continue

            if line.startswith("#") or current_section_key is None:
                continue

            parts = line.split()
            try:
                if current_section_key == "initial_population" and len(parts) >= 3:
                    data["initial_population"].append((int(parts[0]), int(parts[2])))
                elif current_section_key in data["conditions"] and current_state:
                    if current_state in ("before digest", "after digest") and len(parts) >= 2:
                        key = "before" if current_state == "before digest" else "after"
                        data["conditions"][current_section_key][key].append((int(parts[0]), int(parts[1])))
                    elif current_state == "turning points" and parts:
                        n_bp_turn = int(parts[0])
                        data["conditions"][current_section_key]["turning_points"].append(n_bp_turn)
                        if current_section_key == "No Additives":
                            data["turning_points"].append(n_bp_turn)
            except (ValueError, IndexError) as exc:
                print(f"Skipping line: {line!r} ({exc})")

    data["initial_population"].sort()
    data["turning_points"].sort()
    for cond in data["conditions"].values():
        for key in ("before", "after", "turning_points"):
            cond[key].sort()

    return data

def parse_hmfb_data(path: str) -> dict:
    """Parse an HMfB data file with 'theory', 'turning_points', 'before', and 'after' sections."""
    data = {"theory": [], "before": [], "after": [], "turning_points": []}
    current_section = None
    section_keywords = {
        "theory": "theory", "before digest": "before", "before hmfb": "before",
        "turning points": "turning_points", "after digest": "after", "after hmfb": "after"}
    try:
        with open(path, "r", encoding="utf-8") as f:
            for raw in f:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    lower = line.lower()
                    current_section = next((sec for kw, sec in section_keywords.items() if kw in lower), current_section)
                    continue
                if current_section is None:
                    continue
                parts = line.split()
                try:
                    if current_section == "theory" and len(parts) >= 3:
                        data["theory"].append((int(parts[0]), int(parts[2])))
                    elif current_section in ("before", "after") and len(parts) >= 2:
                        data[current_section].append((int(parts[0]), int(parts[1])))
                    elif current_section == "turning_points" and parts:
                        data["turning_points"].append(int(parts[0]))
                except (ValueError, IndexError):
                    print(f"Skipping line in section '{current_section}': {line!r}")
    except FileNotFoundError:
        print(f"Data file not found: {path}")
        raise SystemExit(1)
    except Exception as exc:
        print(f"Could not read file {path}: {exc}")
        raise SystemExit(1)

    for key in data:
        data[key].sort()
    return data

# helper functions
def calc_main_data(n_bp_lim: tuple, turn_points_list: list, min_lk: int) -> tuple:
    """Calculate N_bp, h, and Lk values for sawtooth or theory lines."""
    n_bps, hs, lks = [], [], []
    cur_lk = min_lk
    tp_idx = 0

    while tp_idx < len(turn_points_list) and n_bp_lim[0] > turn_points_list[tp_idx]:
        cur_lk += 1
        tp_idx += 1

    for n_bp in range(n_bp_lim[0], n_bp_lim[1] + 1):
        if tp_idx < len(turn_points_list) and n_bp > turn_points_list[tp_idx]:
            cur_lk += 1
            tp_idx += 1

        h_val = n_bp / cur_lk if cur_lk != 0 else np.nan
        n_bps.append(n_bp)
        hs.append(h_val)
        lks.append(cur_lk)

        if tp_idx < len(turn_points_list) and n_bp == turn_points_list[tp_idx]:
            h_upper = n_bp / (cur_lk + 1) if (cur_lk + 1) != 0 else np.nan
            n_bps.append(n_bp)
            hs.append(h_upper)
            lks.append(cur_lk + 1)

    return n_bps, hs, lks

def get_condition_splits(before_points: list) -> dict:
    """Return {n_bp: (h_min, h_max)} for all N_bp with multiple Lk values."""
    h_by_nbp = defaultdict(list)
    for n_bp, lk in before_points:
        if lk != 0:
            h_by_nbp[n_bp].append(n_bp / lk)
    return {n_bp: (min(vals), max(vals)) for n_bp, vals in h_by_nbp.items() if len(vals) > 1}

def safe_h2phi(h: float) -> float:
    """Convert helical repeat h to twist angle phi (degrees)."""
    arr = np.asarray(h, dtype=float)
    return 360.0 / np.where(arr != 0, arr, 1e-9)

def safe_phi2h(phi: float) -> float:
    """Convert twist angle phi (degrees) to helical repeat h."""
    arr = np.asarray(phi, dtype=float)
    return 360.0 / np.where(arr != 0, arr, 1e-9)

def nbp_to_radius(n_bp: float) -> float:
    """Convert N_bp to radius of curvature (nm)."""
    return np.asarray(n_bp, dtype=float) * rise_per_bp / (2.0 * np.pi)

def radius_to_nbp(radius: float) -> float:
    """Convert radius of curvature (nm) to N_bp."""
    return np.asarray(radius, dtype=float) * (2.0 * np.pi) / rise_per_bp

def ensure_parent_dir(path: str) -> None:
    """Create the parent directory of path if it does not exist."""
    directory = os.path.dirname(path)
    if directory:
        os.makedirs(directory, exist_ok=True)

# plotting helpers
def _setup_secondary_yaxis(ax: mpl.axes.Axes, ylim: tuple) -> None:
    """Add twist-angle secondary y-axis on the rightmost column."""
    ax2 = ax.secondary_yaxis("right", functions=(safe_h2phi, safe_phi2h))
    ax2.set_ylim(360 / ylim[1], 360 / ylim[0])
    ax2.yaxis.set_major_locator(mpl.ticker.MultipleLocator(2))
    ax2.tick_params(axis="y", which="major", length=TICK.major_length)
    ax2.set_ylabel(r"$\phi°$", labelpad=10, rotation=270)

def _setup_secondary_xaxis(ax: mpl.axes.Axes, xlim: tuple) -> None:
    """Add radius-of-curvature secondary x-axis on the top row."""
    ax3 = ax.secondary_xaxis("top", functions=(nbp_to_radius, radius_to_nbp))
    r_min, r_max = nbp_to_radius(xlim[0]), nbp_to_radius(xlim[1])
    ax3.set_xlim(r_min, r_max)
    r_range = r_max - r_min
    thresholds = [5, 10, 20]
    tick_sizes = [0.5, 1.0, 2.0, 5.0]
    major_tick = tick_sizes[np.searchsorted(thresholds, r_range)]
    ax3.xaxis.set_major_locator(mpl.ticker.MultipleLocator(major_tick))
    ax3.tick_params(axis="x", which="major", length=TICK.major_length, width=LW.tick, direction="out")
    ax3.set_xlabel(r"Radius (nm)", labelpad=5)
    ax3.xaxis.set_major_formatter(StrMethodFormatter("{x:.1f}"))
    ax3.xaxis.set_tick_params(labeltop=True)

def _setup_ticks_and_grid(ax: mpl.axes.Axes, xlim: tuple, ylim: tuple,
                           is_first_row: bool, is_last_row: bool, is_first_col: bool, is_last_col: bool) -> None:
    """Configure tick locators, formatters, labels, and background grid."""
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    if is_first_col:
        ax.set_ylabel(r"$h$ (bp/turn)", labelpad=RCPARAMS["axes.labelpad"])
    else:
        ax.yaxis.set_visible(False)
    if is_last_row:
        ax.set_xlabel(r"$N_{\mathregular{bp}}$", labelpad=RCPARAMS["axes.labelpad"])
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_tick_params(which="major", size=TICK.major_length, width=LW.tick, direction="out", top=False)
        axis.set_tick_params(which="minor", size=TICK.minor_length, width=LW.tick, direction="out", top=False)
        axis.set_minor_formatter(NullFormatter())
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, pos: "" if (x == 50 and is_last_col) else f"{x:.0f}"))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(0.5))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.1))
    ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.1f}"))
    ax.tick_params(axis="x", pad=TICK.pad_major)

    if not is_first_row and is_first_col:
        yticks = [t for t in ax.yaxis.get_major_locator().tick_values(*ax.get_ylim())
                  if ax.get_ylim()[0] <= t <= ax.get_ylim()[1]]
        ax.set_yticks(yticks)
        ax.set_yticklabels([f"{t:.1f}" if t < max(yticks) else "" for t in yticks])

    if not is_last_row:
        ax.xaxis.set_ticklabels([])
        ax.tick_params(axis="x", which="both", bottom=False)

    for y_val in (10.5, 11.5, 12.5):
        ax.axhline(y=y_val, linestyle="--", linewidth=LW.grid, color="lightgray", zorder=-20)
    ax.grid(which="both", axis="x", linestyle="--", linewidth=LW.grid, color="lightgray", zorder=-20)

def setup_panel(ax: mpl.axes.Axes, xlim: tuple, ylim: tuple, is_last_row: bool,
                is_first_col: bool, is_last_col: bool, is_first_row: bool) -> None:
    """Configure axes, ticks, secondary axes, labels, and background grid."""
    if is_last_col:
        _setup_secondary_yaxis(ax, ylim)
    if is_first_row:
        _setup_secondary_xaxis(ax, xlim)
    _setup_ticks_and_grid(ax, xlim, ylim, is_first_row, is_last_row, is_first_col, is_last_col)

def plot_panel_background(ax: mpl.axes.Axes, hs_theory_all_lks: list, lk_lim_theory: tuple, n_bp_lim_theory: tuple) -> None:
    """Draw background theory dots."""
    xs = range(n_bp_lim_theory[0], n_bp_lim_theory[1] + 1)
    for i in range(lk_lim_theory[1] - lk_lim_theory[0] + 1):
        ax.plot(xs, hs_theory_all_lks[i], "o", markerfacecolor=colors["theory_dots"],
                markeredgecolor="none", markersize=MS.medium, zorder=2)

def _get_split_colors(condition_label: str) -> tuple:
    """Return (light_color, dark_color) for split lines based on condition."""
    if condition_label == "No Additives":
        return colors["light_green_split"], colors["dark_green_split"]
    if "HMfB" in condition_label:
        return colors.get("hmfb_light_split", "lightcoral"), colors.get("hmfb_band", "red")
    return colors.get("etbr_light_split", "lightblue"), colors.get("etbr_dark_split", "blue")

def plot_reference_splits(ax: mpl.axes.Axes, no_additives_all_splits: dict, panel_type: str) -> None:
    """Draw faint 'No Additives' reference split lines on non-NA panels."""
    custom = custom_splits.get(panel_type, set())
    for n_bp, (h_min, h_max) in no_additives_all_splits.items():
        color = colors["dark_green_split"] if n_bp in custom else colors["light_green_split"]
        ax.plot([n_bp, n_bp], [h_min, h_max], color=color, linewidth=LW.split, alpha=0.7, zorder=1)

def plot_condition_splits(ax: mpl.axes.Axes, condition_label: str, split_points: dict, panel_type: str) -> None:
    """Draw condition-specific split lines (green / blue / red)."""
    light_color, dark_color = _get_split_colors(condition_label)
    exclusion_set = next((n_bp_set for key, n_bp_set in condition_split_exclusions.items() if key in condition_label), set())

    for n_bp, (h_min, h_max) in split_points.items():
        if n_bp in exclusion_set:
            continue
        if condition_label == "No Additives":
            if n_bp in green_splits:
                color = colors["green_split"]
            elif n_bp in custom_splits.get(panel_type, set()):
                color = dark_color
            else:
                color = light_color
        else:
            dark_nbps = next((nbp_list for key, nbp_list in condition_splits.items() if key in condition_label), [])
            color = dark_color if n_bp in dark_nbps else light_color
        ax.plot([n_bp, n_bp], [h_min, h_max], color=color, linewidth=LW.split, alpha=0.7, zorder=3)

def _point_color(n_bp: int, lk: int, condition_label: str, panel_type: str, na_bands: set) -> str:
    """Return the marker color for a data point."""
    is_new   = condition_label != "No Additives" and bool(na_bands) and (n_bp, lk) not in na_bands
    is_faint = (n_bp, lk, condition_label, panel_type) in faint_points
    if is_new:
        return colors["hmfb_band"] if "HMfB" in condition_label else colors["etbr_band"]
    if is_faint:
        return "darkgray"
    return "black"

def plot_condition_data(ax: mpl.axes.Axes, condition_label: str, before_points: list,
                        after_points_set: set, panel_type: str, no_additives_data: dict = None) -> None:
    """Plot primary data circles for a before/after panel."""
    na_bands = set(no_additives_data.get("before", [])) if no_additives_data else set()

    for n_bp, lk in before_points:
        if lk == 0:
            continue
        h_val    = n_bp / lk
        pt_color = _point_color(n_bp, lk, condition_label, panel_type, na_bands)
        in_after = (n_bp, lk) in after_points_set
        keep     = panel_type == "before" or in_after
        ax.plot(n_bp, h_val, "o",
                markersize=MS.large if keep else MS.hollow,
                markerfacecolor=pt_color if keep else "white",
                markeredgecolor="none"   if keep else pt_color,
                markeredgewidth=0        if keep else LW.edge,
                zorder=10)

def _add_before_arrows(ax: mpl.axes.Axes, condition_label: str) -> None:
    """Add annotation arrows to the 'before' panel for specific conditions."""
    for key, arrows in before_arrows.items():
        if key in condition_label:
            for xy, xytext, color, lw in arrows:
                ax.annotate("", xy=xy, xytext=xytext, xycoords="data",
                            arrowprops=dict(arrowstyle="->", color=color, lw=lw, alpha=0.7))

def create_combo_legend(fig: mpl.figure.Figure) -> None:
    """Build legend handles (legend display currently disabled)."""
    handles = [
        Line2D([], [], linestyle="None", marker="o", markersize=MS.large, markerfacecolor="black", markeredgecolor="none", label="Observed circle"),
        Line2D([], [], linestyle="None", marker="o", markersize=MS.hollow, markerfacecolor="white", markeredgecolor="black", markeredgewidth=LW.edge, label="Digested"),
        Line2D([], [], linestyle="None", marker="o", markersize=MS.large, markerfacecolor="darkgray", markeredgecolor="darkgray", label="Faint observed circle"),
        Line2D([], [], linestyle="None", marker="o", markersize=MS.large, markerfacecolor=colors["etbr_band"], markeredgecolor=colors["etbr_band"], label="Only with EtBr"),
        Line2D([], [], linestyle="None", marker="o", markersize=MS.large, markerfacecolor=colors["hmfb_band"], markeredgecolor=colors["hmfb_band"], label="Only with HMfB")
        ]
    # Uncomment to enable:
    # fig.legend(handles=handles, loc="lower center", ncol=5,
    #            frameon=False, bbox_to_anchor=(0.5, -0.15))


# figure generator
def create_figure(condition_list: list, title: str, output_file: str, state: dict) -> None:
    """Create and save a 2-panel (Before/After) figure for each condition row."""
    n_rows = len(condition_list)
    if n_rows == 0:
        print(f"No conditions for {output_file}, skipping.")
        return

    print(f"Creating {output_file} with {n_rows} row(s)...")
    fig, axes = plt.subplots(n_rows, 2, figsize=fig_size_inches(), sharex=False, sharey=False,
                             gridspec_kw={"wspace": 0.01, "hspace": 0.03})
    if n_rows == 1:
        axes = axes.reshape(1, -1)

    fig.suptitle(title, fontsize=RCPARAMS["font.size"], y=1.0 - (0.05 / n_rows))
    no_additives_data = state["conditions"].get("No Additives")

    for row, condition_label in enumerate(condition_list):
        if condition_label not in state["conditions"]:
            print(f"Condition '{condition_label}' not found. Skipping row.")
            continue

        print(f"Row {row + 1}/{n_rows}: {condition_label}")
        ax_before, ax_after = axes[row, 0], axes[row, 1]
        is_first_row = row == 0
        is_last_row  = row == n_rows - 1

        cond_data      = state["conditions"][condition_label]
        before_points  = cond_data["before"]
        after_points_s = set(cond_data["after"])
        split_points   = get_condition_splits(before_points)
        row_ylim       = state["custom_ylims"].get(condition_label, state["ylim"])

        for ax, panel_type in ((ax_before, "before"), (ax_after, "after")):
            is_before = panel_type == "before"
            setup_panel(ax, state["xlim"], row_ylim, is_last_row,
                        is_first_col=is_before, is_last_col=not is_before,
                        is_first_row=is_first_row)
            plot_panel_background(ax, state["hs_theory_all_lks"], state["lk_lim_theory"], state["xlim"])
            if condition_label != "No Additives":
                plot_reference_splits(ax, state["no_additives_all_splits"], panel_type)
            plot_condition_splits(ax, condition_label, split_points, panel_type)
            plot_condition_data(ax, condition_label, before_points, after_points_s, panel_type, no_additives_data)

        _add_before_arrows(ax_before, condition_label)
        ax_before.text(-0.25, 0.5, condition_label, transform=ax_before.transAxes,
                       ha="center", va="center", rotation=90, fontsize=RCPARAMS["font.size"],
                       weight="bold", color=condition_label_colors.get(condition_label, "black"))

    print(f"Saving to {output_file}...")
    ensure_parent_dir(output_file)
    plt.savefig(output_file, dpi=600, bbox_inches="tight", pad_inches=0, facecolor="white")
    print(f"Saved {output_file}.")
    plt.close(fig)

# main
def main() -> None:
    parser = argparse.ArgumentParser(
        prog="plot_combo_paper.py",
        description="Generates combo sawtooth plots for BAL-31, EtBr, and HMFB.",
    )
    parser.add_argument("output", help="Output path prefix (e.g. '260403_NA').")
    parser.add_argument("-x", "--xlim", dest="xlim", help="Comma-separated x-limits (e.g. '50,106').")
    parser.add_argument("-y", "--ylim", dest="ylim", help=f"Comma-separated y-limits. Default: {default_y_lim}.")
    parser.add_argument("--bal31-data", dest="bal31_data", required=True,
                        help="BAL-31/EtBr data file (must contain 'No Additives').")
    parser.add_argument("--hmfb-data", dest="hmfb_data", required=True, help="HMfB data file.")
    args = parser.parse_args()

    for path, label in ((args.bal31_data, "BAL-31"), (args.hmfb_data, "HMFB")):
        if not os.path.isfile(path):
            print(f"{label} data file not found: {path}")
            raise SystemExit(1)

    print(f"Parsing data: {args.bal31_data}")
    bal31_data = parse_bal31_data(args.bal31_data)
    conditions = bal31_data["conditions"]
    na_name        = next((k for k in conditions if "no additive" in k.lower()), None)
    etbr_12_5_name = next((k for k in conditions if "12.5" in k), None)
    if na_name is None:
        print("'No Additives' condition not found. Cannot proceed.")
        raise SystemExit(1)
    if etbr_12_5_name is None:
        print("'EtBr 12.5 ug/mL' condition not found. Cannot proceed.")
        raise SystemExit(1)
    print(f"Found: {na_name!r}, {etbr_12_5_name!r}")

    hmfb_name = "HMfB"
    print(f"Parsing HMFB data: {args.hmfb_data}")
    conditions[hmfb_name] = parse_hmfb_data(args.hmfb_data)
    print(f"Added '{hmfb_name}' to conditions.")

    initial_points = bal31_data["initial_population"]
    turn_points    = bal31_data["turning_points"]
    min_lk   = min(lk  for _, lk  in initial_points)
    min_n_bp = min(nbp for nbp, _ in initial_points)
    max_n_bp = max(nbp for nbp, _ in initial_points)

    xlim = tuple(int(v)   for v in args.xlim.split(",")) if args.xlim else (min_n_bp - 5, max_n_bp + 5)
    ylim = tuple(float(v) for v in args.ylim.split(",")) if args.ylim else default_y_lim
    print(f"X-limits: {xlim}  Y-limits: {ylim}")

    lk_lim_theory = (int(np.floor(xlim[0] / h_canonical)), int(np.ceil(xlim[1] / h_canonical)))

    print("Calculating theory dots...")
    hs_theory_all_lks = [calc_main_data(xlim, [], lk)[1]
                         for lk in range(lk_lim_theory[0], lk_lim_theory[1] + 1)]

    print("Calculating polynomial fit...")
    fit_x_pts, fit_y_pts = [], []
    if len(turn_points) > 5:
        lk_0s  = [min_lk + i + 0.5 for i in range(len(turn_points))]
        h_0s   = [nbp / lk0 for nbp, lk0 in zip(turn_points, lk_0s)]
        coeffs = np.polyfit(turn_points, h_0s, 5)
        fit_x_pts = np.linspace(xlim[0], xlim[1], 400)
        fit_y_pts = np.poly1d(coeffs)(fit_x_pts)
        print("Fit done.")
    else:
        print(f"Only {len(turn_points)} turning points; need > 5 for fit.")

    print("Calculating No Additives reference splits...")
    no_additives_all_splits = get_condition_splits(conditions[na_name]["before"])
    print(f"Found {len(no_additives_all_splits)} reference splits.")

    root, ext = os.path.splitext(args.output)
    if not ext:
        ext = ".svg"

    state = {
        "conditions": conditions,
        "xlim": xlim,
        "ylim": ylim,
        "custom_ylims": {
            na_name:        (9.9, 13.5),
            etbr_12_5_name: (9.9, 14.0),
            hmfb_name:      (9.9, 13.5),
        },
        "hs_theory_all_lks":       hs_theory_all_lks,
        "lk_lim_theory":           lk_lim_theory,
        "min_lk":                  min_lk,
        "fit_x_pts":               fit_x_pts,
        "fit_y_pts":               fit_y_pts,
        "no_additives_all_splits": no_additives_all_splits,
    }

    plots = [
        (f"{root}_etbr_hmfb_combo",          [etbr_12_5_name, hmfb_name]),
        (f"{root}_no_additives_before_after", [na_name]),
    ]
    for i, (base, cond_list) in enumerate(plots, 1):
        print(f"\nGenerating Plot {i}...")
        for fmt in (".png", ".svg"):
            create_figure(condition_list=cond_list, title="", output_file=f"{base}{fmt}", state=state)
    print("\nAll combo plots generated.")

if __name__ == "__main__":
    main()
