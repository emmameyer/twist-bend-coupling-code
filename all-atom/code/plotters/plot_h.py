"""
plot_h.py

plot helical repeat (h) over time for one or more simulations (simple minicircular DNA or minicircular DNA with proteins)
structured to support either case, two types of input (i have it set up to take my twist xvg file from the process_twist.py script):

    1. simple input: a single xvg file per simulation (no proteins). Use the following flag:
        --sim path/to/twist_internal_avg.xvg --label "No Histone (for example)"
    
    2. split input: pair of xvg files per simulation (one for the region of DNA interaction with the protein and one for the non-interaction region).
        Use the following flags:
        --sim path/to/interacting/twist_internal_avg.xvg \
                path/to/noninteracting/twist_internal_avg.xvg \
        --label "Nick + Histone (for example)"

The --sim and --label flags can be repeated for as many simulations as needed.

Example usage:

python3 plot_h.py \
--sim data/nick_no_histone/twist_internal_avg.xvg \
--label "Nick, No Histone" \
--sim data/interacting/twist_internal_avg.xvg \
        data/noninteracting/twist_internal_avg.xvg \
--label "Nick with Histone" \
--sim data/no_nick_no_histone/twist_internal_avg.xvg \
--label "No Nick, No Histone" \
--plot_name whatever_you_call_your_plot \
--window 50


python3 plot_h.py \
--sim /Users/emma/Documents/kent/5t5k/paper_specific/nick_no_histone/twist_internal_avg.xvg \
--label "No Histone" \
--sim /Users/emma/Documents/kent/5t5k/paper_specific/nick_histone/interacting_noninteracting_separate/interacting/twist_internal_avg.xvg \
/Users/emma/Documents/kent/5t5k/paper_specific/nick_histone/interacting_noninteracting_separate/noninteracting/twist_internal_avg.xvg \
--label "Histone" \
--plot_name h_overtime_histone_comparison \
--window 50

This script does NOT auto detect interacting/noninteracting regions, you have to use two inputs. The input files should already separate interacting and noninteracting.

"""

import argparse
import sys
import numpy as np
import matplotlib.pyplot as plt
from pyparsing import line
# i have these commented out, we liked the rolling average better. uncomment if you want data fit with curve_fit or linregress
# from scipy.optimize import curve_fit
# from scipy.stats import linregress

from fig_style import apply_style, LW, TICK
apply_style()

# aesthetics

colors = [
    "blue", "green", "red", "orange", "purple",
    "brown", "pink", "gray", "olive", "cyan",
]

# option to change linestyle for the different cases
linestyles = {
    "simple":         "-",   # single file
    "interacting":    "-",   # histone interacting region
    "noninteracting": "-",   # noninteracting (free) region
}


# arg parsing

class SimAction(argparse.Action):
    """Collect each --sim (simulation) ... into a list of lists."""
    def __call__(self, parser, namespace, values, option_string=None):
        items = getattr(namespace, self.dest, None) or []
        items.append(values)
        setattr(namespace, self.dest, items)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot helical repeat over time; supports histone interacting simulations.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--sim", dest="sims", nargs="+", action=SimAction, metavar="XVG",
        required=True,
        help="One xvg file (simple) or two xvg files (interacting / noninteracting) per simulation being analyzed. "
             "Repeat --sim for each simulation.",
    )
    parser.add_argument(
        "--label", dest="labels", action="append", metavar="LABEL",
        default=None,
        help="Display label for the following --sim. "
             "Must be given once per --sim (in the same order).",
    )
    parser.add_argument(
        "--plot_name", dest="plot_name", default="output", metavar="NAME",
        help="File name for the output plot, default: 'output'.",
    )
    parser.add_argument(
        "--window", dest="window", type=int, default=50, metavar="N",
        help="Rolling average window size in frames (default: 50).",
    )
    args = parser.parse_args()

    n = len(args.sims)                          # total number of simulations

    # make sure only 1 or 2 files are given for each simulation and that the number of labels matches as well
    for i, paths in enumerate(args.sims):
        if len(paths) not in (1, 2):
            parser.error(
                f"--sim entry {i+1} has {len(paths)} file(s); "
                "expected 1 (simple) or 2 (interacting / noninteracting)."
            )

    # auto generates a label (Sim #)
    if args.labels is None:
        args.labels = [f"Sim {i+1}" for i in range(n)]
    elif len(args.labels) != n:
        parser.error(
            f"{len(args.labels)} --label(s) given but {n} --sim group(s) found; "
            "provide exactly one --label per --sim."
        )

    return args


# load data

def load_xvg(path):
    """Return (time_ns, helical_repeat_bp_per_turn)."""
    try:
        data = np.loadtxt(path, comments=["@", "#"])
    except Exception as e:
        sys.exit(f"Error reading {path}: {e}")
    time = data[:, 0] / 1000          # ps to ns
    h    = 360.0 / data[:, 1]         # deg to bp/turn
    return time, h


# rolling average

def rolling_average(time, h, window):
    """Return (time, rolling_avg) using pandas rolling mean with min_periods=1.
    min_periods=1 means the average is computed from whatever data is available
    at the edges, so the output always covers the full time range with no trimming (looks better on edges)."""
    import pandas as pd
    h_avg = pd.Series(h).rolling(window=window, center=True, min_periods=1).mean().to_numpy()
    return time, h_avg


# uncomment for fitting instead of rolling avg
#
# def detect_fit(time, h):
#     slope, intercept, _, p, _ = linregress(time, h)
#     initial = np.mean(h[: max(1, len(h) // 10)])
#     final   = np.mean(h[-max(1, len(h) // 10):])
#     delta   = final - initial
#     if abs(delta) > 0.1:
#         fit_type = "exponential"
#     elif p < 0.05:
#         fit_type = "linear"
#     else:
#         fit_type = "equilibrated"
#     return fit_type, slope, intercept, p, delta
#
#
# def add_fit(ax, time, h, fit_type, slope, intercept, color, linestyle, label):
#     """Draw the fit line and return a legend-ready label string."""
#     t_dense = np.linspace(time[0], time[-1], 2000)
#
#     if fit_type == "exponential":
#         def exp_decay(t, A, tau, C):
#             return A * np.exp(-t / tau) + C
#         try:
#             p0 = (h[0] - h[-1], (time[-1] - time[0]) / 5, h[-1])
#             popt, _ = curve_fit(exp_decay, time, h, p0=p0, maxfev=10_000)
#             _, tau_fit, _ = popt
#             print(f"    Exp fit: τ = {tau_fit:.2f} ns")
#             ax.plot(t_dense, exp_decay(t_dense, *popt),
#                     color=color, linewidth=LW.fit, linestyle=linestyle,
#                     label=f"{label} (τ = {tau_fit:.1f} ns)")
#         except Exception as e:
#             print(f"    Exponential fit failed: {e}")
#
#     elif fit_type == "linear":
#         ax.plot(t_dense, slope * t_dense + intercept,
#                 color=color, linewidth=LW.fit, linestyle=linestyle,
#                 label=f"{label} (slope={slope:.2e})")
#
#     else:  # equilibrated: mean ± std
#         mean_val = np.mean(h)
#         std_val  = np.std(h)
#         ax.axhline(mean_val, color=color, linewidth=LW.fit, linestyle=linestyle,
#                    label=f"{label} ({mean_val:.3f}±{std_val:.3f} bp/turn)")


# main

def main():
    args = parse_args()

    fig, ax = plt.subplots(figsize=(2.5, 2.0), tight_layout=True)
    all_time, all_h = [], []

    color_idx = 0  # increment for every plotted series (just goes thru the list of colors defined above in order)

    for sim_idx, (paths, label) in enumerate(zip(args.sims, args.labels)):
        is_split = len(paths) == 2

        if is_split:
            # if there are two files, the first is hardcoded as interacting and the second as noninteracting
            interacting_path,    = paths[:1]
            noninteracting_path, = paths[1:]

            for path, kind in [
                (interacting_path,    "interacting"),
                (noninteracting_path, "noninteracting"),
            ]:
                color = colors[color_idx % len(colors)]
                color_idx += 1

                time, h = load_xvg(path)
                all_time.extend(time)
                all_h.extend(h)

                sublabel = f"{'Interacting with' if kind == 'interacting' else 'Not interacting with'} {label}"

                # raw data (transparent)
                ax.plot(time, h, color=color, linewidth=LW.data,
                        alpha=0.35, linestyle=linestyles[kind])

                # rolling average
                time_avg, h_avg = rolling_average(time, h, args.window)
                ax.plot(time_avg, h_avg, color=color, linewidth=LW.fit,
                        linestyle=linestyles[kind], label=sublabel)

                # uncomment for fitting instead of rolling avg
                # fit_type, slope, intercept, p_val, delta = detect_fit(time, h)
                # print(f"[{sublabel}] {fit_type}  Δh={delta:.3f} bp/turn  p={p_val:.3f}")
                # add_fit(ax, time, h, fit_type, slope, intercept, color, ls, sublabel)

        else:
            color = colors[color_idx % len(colors)]
            color_idx += 1
            kind = "simple"

            # simple case (no histone / don't care about interacting vs noninteracting)
            path = paths[0]
            time, h = load_xvg(path)
            all_time.extend(time)
            all_h.extend(h)

            # raw data (transparent)
            ax.plot(time, h, color=color, linewidth=LW.data,
                    alpha=0.35, linestyle=linestyles[kind])

            # rolling average
            time_avg, h_avg = rolling_average(time, h, args.window)
            ax.plot(time_avg, h_avg, color=color, linewidth=LW.fit,
                    linestyle=linestyles[kind], label=label)

            # uncomment for fitting instead of rolling avg
            # fit_type, slope, intercept, p_val, delta = detect_fit(time, h)
            # print(f"[{label}] {fit_type}  Δh={delta:.3f} bp/turn  p={p_val:.3f}")
            # add_fit(ax, time, h, fit_type, slope, intercept, color, ls, label)

    # B-DNA reference
    #ax.axhline(10.5, color="black", linewidth=LW.grid, linestyle="--",
    #           label="B-DNA (10.5 bp/turn)")

    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("h (bp/turn)")
    ax.set_xlim(0, max(all_time))
    ax.set_ylim(9.7, 11.8)
    ax.set_xticks(np.arange(0, max(all_time) + 1, 5))
    #ax.legend(fontsize=6)   # uncomment if you want the legend

    out = "{}.svg".format(args.plot_name)
    fig.savefig(out, bbox_inches="tight", dpi=600, facecolor="white", pad_inches=0.02)
    plt.close(fig)
    print(f"\nSaved → {out}")


if __name__ == "__main__":
    main()
