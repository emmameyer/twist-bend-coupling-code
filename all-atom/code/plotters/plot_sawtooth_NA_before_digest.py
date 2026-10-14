"""
plot_sawtooth_NA_before_digest.py

plots just the before panel of the no additives Soumya gel data, with a polynomial fit
to use:
python3 plot_sawtooth_NA_before_digest.py 03_27_2026_na_before.svg \
--data-file bal31_data.txt \
--xlim 50,106 \
--ylim 9.0,13.5
"""

import argparse
import os
import re
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import defaultdict
from matplotlib.ticker import NullFormatter, StrMethodFormatter
from typing import Dict
from fig_style import RCPARAMS, MS, LW, fig_size_inches

# constants
h_canonical  = 10.5
rise_per_bp  = 0.34
default_ylim = (9.7, 14.5)
colors = {
    "fit_line":    "green",
    "green_split": "#035803",
    "light_green": "green",
    "theory_dots": "lightgray",
    "faint":       "darkgray",
}
green_splits  = {73, 83, 92, 93, 103}
custom_splits = {63, 83, 93}
split_midpoints = [
    (53, 11.7777777777777),
    (63, 11.4545454545),
    (73, 11.2307692308),
    (83, 11.066666667),
    (93, 10.9411764706),
    (103, 10.8421052632),
]

# fig size
fig_width_mm  = 53
fig_height_mm = 37

# secondary axis conversions
def h2phi(h):
    h = np.asarray(h, dtype=float)
    return 360.0 / np.where(h != 0, h, 1e-9)

def phi2h(phi):
    phi = np.asarray(phi, dtype=float)
    return 360.0 / np.where(phi != 0, phi, 1e-9)

def nbp_to_radius(n_bp):
    return np.asarray(n_bp, dtype=float) * rise_per_bp / (2.0 * np.pi)

def radius_to_nbp(r):
    return np.asarray(r, dtype=float) * (2.0 * np.pi) / rise_per_bp

# data parsing
def parse_digest_data(path: str) -> Dict:
    data = {"initial_population": [], "turning_points": [], "conditions": {}}
    section_key   = None
    current_state = None

    with open(path, encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.lower().startswith("# theory"):
                continue
            if line.startswith("# N_bp") and "h_expected" in line:
                section_key   = "initial_population"
                current_state = None
                continue
            m = re.match(r"# (.*?) (before digest|after digest|turning points)",
                         line, re.IGNORECASE)
            if m:
                label = m.group(1).strip()
                current_state = m.group(2).strip().lower()
                if label.lower() in ("no additives", "no additive"):
                    label = "No Additives"
                if label not in data["conditions"]:
                    data["conditions"][label] = {"before": [], "after": [], "turning_points": []}
                section_key = label
                continue
            if line.startswith("#") or section_key is None:
                continue
            parts = line.split()
            try:
                if section_key == "initial_population" and len(parts) >= 3:
                    data["initial_population"].append((int(parts[0]), int(parts[2])))
                elif section_key in data["conditions"] and current_state:
                    if current_state in ("before digest", "after digest") and len(parts) >= 2:
                        key = "before" if current_state == "before digest" else "after"
                        data["conditions"][section_key][key].append(
                            (int(parts[0]), int(parts[1]))
                        )
                    elif current_state == "turning points" and len(parts) >= 1:
                        tp = int(parts[0])
                        data["conditions"][section_key]["turning_points"].append(tp)
                        if section_key == "No Additives":
                            data["turning_points"].append(tp)
            except (ValueError, IndexError):
                pass
    data["initial_population"].sort()
    data["turning_points"].sort()
    for c in data["conditions"].values():
        c["before"].sort()
        c["after"].sort()
        c["turning_points"].sort()

    return data

# find splits (multiple Lk values for same nbp)
def get_splits(before_points):
    h_by_nbp = defaultdict(list)
    for n_bp, lk in before_points:
        if lk:
            h_by_nbp[n_bp].append(n_bp / lk)
    return {n: (min(hs), max(hs)) for n, hs in h_by_nbp.items() if len(hs) > 1}

# setup
def setup_panel(ax, xlim, ylim):
    # Right y-axis: uncomment if you want to show phi values
    #ax2 = ax.secondary_yaxis("right", functions=(h2phi, phi2h))
    #ax2.set_ylim(360 / ylim[1], 360 / ylim[0])
    #ax2.yaxis.set_major_locator(mpl.ticker.MultipleLocator(2))
    #ax2.tick_params(axis="y", which="major", length=2, width=0.4, pad=1.5)
    #ax2.set_ylabel(r"$\phi°$", labelpad=5, rotation=270, fontsize=6)

    # Top x-axis: radius (nm)
    ax3 = ax.secondary_xaxis("top", functions=(nbp_to_radius, radius_to_nbp))
    r_min, r_max = nbp_to_radius(xlim[0]), nbp_to_radius(xlim[1])
    ax3.set_xlim(r_min, r_max)
    r_range = r_max - r_min
    major_tick = 0.5 if r_range < 5 else (1.0 if r_range < 10 else (2.0 if r_range < 20 else 5.0))
    ax3.xaxis.set_major_locator(mpl.ticker.MultipleLocator(major_tick))
    ax3.set_xlabel(r"Radius (nm)", labelpad=2, fontsize=6)
    ax3.xaxis.set_major_formatter(StrMethodFormatter("{x:.1f}"))
    ax3.xaxis.set_tick_params(labeltop=True)
    ax.set_ylabel(r"$h$ (bp/turn)", labelpad=2)
    ax.set_xlabel(r"$N_{\mathregular{bp}}$", labelpad=2)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    for axis in [ax.xaxis, ax.yaxis]:
        axis.set_minor_formatter(NullFormatter())
    ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))
    ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
    ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(1.0))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(0.5))
    ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.1f}"))

    #ax.axhline(h_canonical, linestyle="-", linewidth=1.0, color="green", zorder=2)

    for y_val in [10.5, 11.5, 12.5]:
        ax.axhline(y_val, linestyle="--", linewidth=LW.grid, color="lightgray", zorder=0, alpha=0.7)
    ax.grid(which="both", axis="x", linestyle="--", linewidth=LW.grid, color="lightgray", zorder=0, alpha=0.7)

# plotting functions
def plot_theory_background(ax, xlim, lk_lim):
    ns = list(range(xlim[0], xlim[1] + 1))
    for lk in range(lk_lim[0], lk_lim[1] + 1):
        hs = [n / lk if lk else np.nan for n in ns]
        ax.plot(ns, hs, "o",
                markersize=MS.medium, markerfacecolor=colors["theory_dots"],
                markeredgecolor="none", zorder=1)

def plot_splits(ax, split_points):
    for n_bp, (h_min, h_max) in split_points.items():
        color = colors["green_split"] if (n_bp in green_splits or n_bp in custom_splits) \
                else colors["light_green"]
        ax.plot([n_bp, n_bp], [h_min, h_max],
                color=color, linewidth=LW.split, alpha=0.5, zorder=1)

def plot_data_points(ax, before_points):
    for n_bp, lk in before_points:
        if lk == 0:
            continue
        ax.plot(n_bp, n_bp / lk, "o",
                markersize=MS.large, markerfacecolor="black",
                markeredgecolor="none", zorder=10)

def plot_fit(ax, xlim, turn_points, min_lk):
    lk_0s = [min_lk + i + 0.5 for i in range(len(turn_points))]
    h_0s  = [N / L for N, L in zip(turn_points, lk_0s)]
    if len(turn_points) > 5:
        poly = np.poly1d(np.polyfit(turn_points, h_0s, 5))
        xs   = np.linspace(xlim[0], xlim[1], 400)
        ax.plot(xs, poly(xs), "-", color=colors["fit_line"], linewidth=1.0, zorder=2)

def plot_midpoints(ax):
    for n_bp, h_val in split_midpoints:
        ax.plot(n_bp, h_val, "x",
                color="green", markersize=MS.cross, markeredgewidth=LW.tick, zorder=15)

# main
def main():
    parser = argparse.ArgumentParser(description="No Additives Gel, Before Digest sawtooth")
    parser.add_argument("output",      help="Output SVG path")
    parser.add_argument("--data-file", required=True, dest="data_file")
    parser.add_argument("--xlim",      default=None)
    parser.add_argument("--ylim",      default=None)
    args = parser.parse_args()

    if not os.path.isfile(args.data_file):
        print(f"Data file not found: {args.data_file}")
        raise SystemExit(1)

    print(f"Parsing {args.data_file} …")
    digest = parse_digest_data(args.data_file)

    initial_points = digest["initial_population"]
    turn_points    = digest["turning_points"]
    conditions     = digest["conditions"]

    if "No Additives" not in conditions:
        print("'No Additives' condition not found in data file.")
        raise SystemExit(1)

    before_points = conditions["No Additives"]["before"]
    min_lk   = min(lk for _, lk in initial_points)
    min_n_bp = min(n  for n,  _ in initial_points)
    max_n_bp = max(n  for n,  _ in initial_points)

    xlim = ((int(args.xlim.split(",")[0]), int(args.xlim.split(",")[1]))
            if args.xlim else (min_n_bp - 5, max_n_bp + 5))
    ylim = ((float(args.ylim.split(",")[0]), float(args.ylim.split(",")[1]))
        if args.ylim else default_ylim)

    lk_lim = (
        int(np.floor(xlim[0] / h_canonical)),
        int(np.ceil(xlim[1]  / h_canonical)),
    )

    plt.rcParams.update(RCPARAMS)
    # uncomment if you want to specify figure size, else its from fig_style
    #fig_w = fig_width_mm  / 25.4
    #fig_h = fig_height_mm / 25.4
    fig, ax = plt.subplots(1, 1, figsize=fig_size_inches(), dpi=600)
    print(ax.get_position())

    setup_panel(ax, xlim, ylim)
    plot_theory_background(ax, xlim, lk_lim)
    plot_splits(ax, get_splits(before_points))
    plot_data_points(ax, before_points)
    plot_fit(ax, xlim, turn_points, min_lk)
    plot_midpoints(ax)
    plt.tight_layout(pad=0.2)

    out = args.output if args.output.endswith(".svg") else args.output + ".svg"
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    plt.savefig(out, dpi=600)
    print(f"Saved -> {out}")
    plt.close(fig)

if __name__ == "__main__":
    main()
