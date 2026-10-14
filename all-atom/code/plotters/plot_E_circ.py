#Plots the free energies of circularization considering 1) twist energy, 2) bending energy, and 3) fluctuations about a closed circle (following some mathematics of Shimada Yamakawa paper).
#It does this both utilizing the polynomial fit for the natural helical repeat from 'make_sawtooth_fig.py' to the experimental datapoints, and assuming hypothetical no TBC result with natural helical repeat fixed at 10.5 bp/turn.
#While the output plot is quite different from the sawtooth plot, internally making both plots are quite the same, so a lot of code is inherited from 'make_sawtooth_fig.py'.

# Note: the notation used in the code follows the Shimada Yamakawa notation, where instead of Tw (Delta_Tw, etc.), they write N (Delta_N, etc.).
# originally written by Tommy and Dr. Portman, modified by Emma Meyer (added some CLI arguments, datapoints to match sawtooth, etc.)

#For usage, run "python plot_E_circ.py -h"

import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from typing import List, Tuple, Callable
from matplotlib.ticker import NullFormatter, StrMethodFormatter
from fig_style import MS, LW, apply_style, fig_size_inches

# parse CLI input parameters
prog = "plot_E_circ.py"
parser = argparse.ArgumentParser(prog = prog, description='''Attempts to plot the free energy of circularization Delta_F = E_bend + E_Tw + Delta_f = E_circ + Delta_f (see manuscript supplementary note 2).
                                 E_bend and E_twist are simply the bend and twist energy necessary to form a closed loop, while delta_F represents the fluctuations in shape about a closed circle (that results in stabilizing the circular state by around 10-20 k_BT iirc).
                                 This utilizes the output data from `make_sawtooth_fig.py` (currently hardcoded) to plot the free energy as a function of DNA length to both the hypothetical no TBC datapoints (i.e., natural helical repeat fixed at 10.5 bp/turn),
                                 and the proposed landscape for the experimental datapoints (using the polynomial fit data for the natural helical repeat as a function of length).''')
parser.add_argument("plot_dir", type=str, nargs=1, help="Output plot directory (name) (.png/.svg).")
parser.add_argument("-x", "--xlim", type=str, nargs=1, dest="xlim", help="optional: Comma-separated x-limits for the plot (inclusive on both ends). Defaults to (60,105), the range of experimental lengths. Note that increasing these ranges further will only affect the theoretical curve unless more experimental data is provided in the hard code of the program.")
parser.add_argument("-y", "--ylim", type=str, nargs=1, dest="ylim", help="optional: Comma-separated y-limits for the plot (inclusive on both ends). Defaults to (10,60).")
parser.add_argument("--no-theory", action = "store_true", dest="no_theory", help="optional: If set, only the experimental datapoints and curve are plotted.")
parser.add_argument("--only-theory", action="store_true", dest="only_theory", help="optional: If set, only the theoretical datapoints and curve are plotted.")

args = parser.parse_args()
plot_dir: str = args.plot_dir[0]
xlim: Tuple[int,int]
ylim: Tuple[float,float]
if args.xlim:
    xlim = (int(args.xlim[0].split(",")[0]), int(args.xlim[0].split(",")[1]))
else:
    xlim = (60, 105)
if args.ylim:
    ylim = (float(args.ylim[0].split(",")[0]), float(args.ylim[0].split(",")[1]))
else:
    ylim = (10, 60)
no_theory: bool = args.no_theory
only_theory: bool = args.only_theory

# experimental datapoint parameters

N_bp_lim: Tuple[int,int] = (60,105)                                                 # min/max values for plotting
turn_points: List[int] = [63, 73, 83, 93, 103]                                      # the values where a split is observed experimentally (Lk_0 shifts from n to n+1)
min_Lk: int = 5                                                                     # minimum Lk_0 value that is plotted (value of Lk from N_bp_lim[0] to turn_points[0])
rel_N_bps: List[int] = [68, 78, 88, 98]                                             # the values (lengths) whre the minicircles are relaxed

#match sawtooth (experimental) data
special_points: dict = {
    5:  [54],
    6:  [63, 64, 65],
    7:  [73, 74, 75],
    8:  [83, 84, 85, 86, 87],
    9:  [92, 93, 94, 95, 96, 97],
    10: [102, 103, 104, 105]
}
# "hollow" points, bp lengths that were digested in the experimental data
digested_points: dict = {
    5:  [55, 56, 57, 58, 59, 60, 61, 62],
    6:  [66, 67, 68, 69, 70, 71, 72],
    7:  [76, 77, 78, 79, 80, 81, 82],
    8:  [88, 89, 90, 91, 92],
    9:  [98, 99, 100, 101, 102]
}
digested_points_theory: dict = {
    6:  [66, 67, 68, 69, 70, 71, 72, 73],
    7:  [76, 77, 78, 79, 80, 81, 82, 83],
    8:  [88, 89, 90, 91, 92, 93],
    9:  [98, 99, 100, 101, 102, 103]
}

# hypothetical/theoretical plot parameters

N_bp_lim_theory: Tuple[int,int] = xlim                                              # min/max values for plotting theortical data
h_canonical: float = 10.5                                                           # canonical value of relaxed helical pitch for linear dsDNA
Lk_lim_theory: Tuple[int,int] = (int(np.round(N_bp_lim_theory[0]/h_canonical)),
                                 int(np.ceil(N_bp_lim_theory[1]/h_canonical)))      # min/max Lk_0 that will be plotted (determined automatically)
turn_points_theory: List[int] = [int(round(h_canonical * (Lk +1./2.)))
                                 for Lk in range(*Lk_lim_theory)]                   # turning points by which Lk_0 shifts from n to n +1 in theory (determined automatically)

# tertiary parameters

Lp: float = 150                                                                     # bending persistence length (in units of bp)
twist_Lp: float = 100/.34                                                           # twisting persistence length (in units of bp)

# functions for calculating the energies

def calc_main_data(N_bp_lim: Tuple[int,int], turn_points: List[int], min_Lk: int) -> Tuple[List[int], List[float], List[float], List[int]]:
    """
    Determines the helical repeat h = N_bp/Lk of duplexes of all integer lengths in provided limits of `N_bp_lim`, 
    given a starting value of `min_Lk` for the lower limit `N_bp_lim[0]` and lengths `turn_points` where the linking number increments so as to minimize the free energy. 
    This creates a "sawtooth" shape, where a change in linking number defines the vertical lines.
    Note 1: If a duplex has length equal to `turn_points[i]` for any `i`, this function computes both the helical pitch at nearest Lk *and* Lk + 1. This can be modified by changing the if statement logic.
    Note 2: If instead you wish to compute the h = N_bp/Lk of all lengths in limits `N_bp_lim` for a nonchanging linking number `min_Lk`, set `turn_points = []`.
    """
    N_bps: List[int] = []; hs: List[float] = []; Lks: List[int] = []                # will hold the three quantities for h = N_bp / Lk
    cur_Lk: int = min_Lk                                                            # starting linking number of 'N_bp_lim[0]' (incremented as needed such that the free energy is approx. minimized)
    for N_bp in range(N_bp_lim[0], N_bp_lim[1] +1):
        cur_index: int = cur_Lk -min_Lk
        if cur_index == len(turn_points) or N_bp <= turn_points[cur_index]:         # record N_bp/cur_Lk
            N_bps.append(N_bp); hs.append(N_bp/cur_Lk); Lks.append(cur_Lk)
        if cur_index < len(turn_points) and N_bp == turn_points[cur_index]:         # if a turning point is reached, increment cur_Lk and record upper values as well
            cur_Lk += 1
            N_bps.append(N_bp); hs.append(N_bp/cur_Lk); Lks.append(cur_Lk)
        
    return N_bps, hs, Lks

def calc_E_Tw(N_bp: int, twist_Lp: float, h_0: float, Lk: float) -> float:
    """
    Twist energy (in units of k_BT) necessary to form a ring (to put the ends in phase)
    E_Tw/k_BT = 2*Pi^2*l_t*Delta_Tw^2/L where l_t is twist persistence length, Delta_Tw is excess twist (from manuscript)
    """
    N_0: float = float(N_bp)/h_0                                                    # equilibrium number of helical turns in duplex length L base pairs 
    Delta_N: float = N_0 - Lk                                                       # excess twist number effectively
    E_Tw: float = 2 * np.pi**2 * twist_Lp / N_bp * abs(Delta_N)**2                  # corresponding twist energy in units of k_BT

    return E_Tw

def calc_E_bend(N_bp: int, Lp: float) -> float:
    """
    Bending energy (in units of k_BT) necessary to form a ring (to bring the ends into proximity) 
    (Thorsten cited Garcia 2007)
    E_bend/k_BT = L*lp/2*kappa^2 where L is duplex length, lp is persistence length, kappa = 2*Pi/L is curvature (from manuscript)
    kappa = 1/r
    """
    E_bend: float = 2 * np.pi**2 * Lp / N_bp

    return E_bend

def calc_Delta_f(N_bp: int, Lp: float, twist_Lp: float, Lk: int, h_0: float) -> float:
    """
    Stabilization energy (in units of k_BT)via small spatial fluctuations about a perfect circle.
    Modified fn. by Dr. Portman 
    """    
    # constants necessary to the calculation valid for SY low L limit
    a_0j: List[float] = [2.784, 2.113, .6558, 1.719, -2.478, 2.588, -1.210, 0.2437]
    a_1j_0: List[float] = [.2639, .1399, -.1131, .6500, -1.1223, 1.0320, -.4601, .0829]
    a_1j_1: List[float] = [-.0383, -.0827, 0.0125, -.2170, .3961, -.3991, .1899, -.0367]

    # compute excess twist from equilibrium to form Lk
    N_0: float = float(N_bp)/h_0
    Delta_N: float = N_0 - Lk

    # compute Poisson's ratio and Kuhn length, convert N_bp to reduced length scale in paper
    sigma: float = Lp/twist_Lp - 1                                                  # Poisson's ratio
    lambda_inv: float = Lp * 2                                                      # Kuhn length (in units of bp)
    red_L: float = N_bp/(lambda_inv)                                                # length of duplex in multiples of Kuhn length (twice persistence length)

    # compute C_0 from eqn. 38 in SY paper
    prefac_C_0: float = 1/(np.sqrt(1 + sigma))
    C_0: float = 0 
    for j in range(len(a_0j)):
        C_0 += prefac_C_0 * a_0j[j] * (Delta_N/(1 + sigma))**(2*j)

    # compute C_1 from eqn. 39 in SY paper
    prefac_C_1: float = 1/(1 + sigma)
    C_1: float = 0
    for j in range(len(a_1j_0)):
        C_1 += prefac_C_1 * (a_1j_0[j] + a_1j_1[j]) * (Delta_N/(1 + sigma))**(2*j)

    # if you use the functional form of G in SY paper, you can find that you can write -ln G in the following form (separated into terms involving E_bend and E_Tw, and the remaining terms encompassing Delta_f/k_BT)
    # Delta_F = -ln G = (E_bend + E_twist)/k+BT + Delta_f/k_BT, where first term is the first term of the exponential in G
    Delta_f: float = -(C_1 + 0.25) * red_L - np.log(C_0) + 13./2.*np.log(red_L)        # fluctuation energy (in units of k_BT)

    return Delta_f

def calc_E_tots(N_bps: List[int], Lks: List[int], fit_fn: Callable[[float],float], Lp: float, twist_Lp: float) -> List[float]:
    """
    Given parallel lists `N_bps` and `Lks` (i.e., from func `calc_main_data`), and a function `fit_fn` for h_0 as a function of `N_bp` (if non-changing, define a function that returns a constant value),
    calculates the total energy for each `N_bp` in `N_bps` E_tot = E_bend + E_Tw + Delta_f = E_circ + Delta_f.
    """
    E_tots: List[float] = []
    for N_bp, Lk in zip(N_bps, Lks):
        h_0: float = fit_fn(N_bp)
        E_Tw: float = calc_E_Tw(N_bp, twist_Lp, h_0, Lk)
        E_bend: float = calc_E_bend(N_bp, Lp)
        Delta_f: float = calc_Delta_f(N_bp, Lp, twist_Lp, Lk, h_0)
        E_tot: float = E_Tw + E_bend + Delta_f
        
        E_tots.append(E_tot)
    
    return E_tots
    

def calc_E_tots_per_Lk(N_bp_lim: Tuple[int, int], unique_Lks: List[int], fit_fn: Callable[[float], float], Lp: float, twist_Lp: float) -> dict:
    """
    For each unique Lk value in `unique_Lks`, computes the total energy E_tot = E_bend + E_Tw + Delta_f
    for every integer N_bp in `N_bp_lim`, holding Lk fixed. Returns a dict mapping each Lk to a tuple
    (N_bps, E_tots) covering the full N_bp range — used to draw the dashed parabola extensions.
    Delta_f is evaluated at Delta_N=0 for each N_bp (i.e., Lk_relaxed = N_bp/h_0) so that the SY polynomial
    remains in its valid range and the dashed lines connect smoothly to the solid curve at junction points.

    This is to extend the parabolas for each Lk.
    """
    result: dict = {}
    for Lk in unique_Lks:
        N_bps_full: List[int] = list(range(N_bp_lim[0], N_bp_lim[1] + 1))
        E_tots_full: List[float] = []
        for N_bp in N_bps_full:
            h_0: float = fit_fn(N_bp)
            Lk_relaxed: float = float(N_bp)/h_0                                 #relaxed Lk (Delta_N=0) used to keep Delta_f in valid SY range
            E_Tw: float = calc_E_Tw(N_bp, twist_Lp, h_0, Lk)
            E_bend: float = calc_E_bend(N_bp, Lp)
            Delta_f: float = calc_Delta_f(N_bp, Lp, twist_Lp, Lk_relaxed, h_0) #Delta_f fixed at Delta_N=0 baseline
            E_tots_full.append(E_Tw + E_bend + Delta_f)
        result[Lk] = (N_bps_full, E_tots_full)
    return result

# determine helical pitches at experimentally lengths 'rel_N_bps' (lightgray circle markers for twist bend coupling experimental h_0 fit)
rel_hs: List[float] = []
rel_Lks: List[int] = []

# determine helical pitches at unrelaxed lengths turn_points (Lk_0 is n + 1/2)
turn_points_Lk_0s: List[float] = [min_Lk + i + 1./2. for i in range(len(turn_points))]
turn_points_h_0s: List[float] = [N_bp/Lk_0 for N_bp, Lk_0 in zip (turn_points, turn_points_Lk_0s)]

# determine the solid N_bp range for each Lk from the main sawtooth data as a (min,max) interval
solid_range_per_Lk: dict = {}
solid_range_per_Lk_theory: dict = {}

# calculate main data, concatenate list of relaxed points, sort some lists, calculate energies, parabolas, 
# separated by only theory/no theory for two different outputs
if not only_theory:
    N_bps: List[int]; hs: List[float]; Lks: List[int]
    N_bps, hs, Lks = calc_main_data(N_bp_lim, turn_points, min_Lk)
    for rel_N_bp in rel_N_bps:
        for i, N_bp in enumerate(N_bps):
            if rel_N_bp == N_bp:
                rel_hs.append(hs[i])
                rel_Lks.append(Lks[i])
    concat_N_bps: List[int] = rel_N_bps + turn_points
    concat_h_0s: List[float] = rel_hs + turn_points_h_0s
    #Sort the lists for plot-related reasons and printing out diagnostics (see: https://stackoverflow.com/questions/5284183/python-sort-list-with-parallel-list)
    zipped_N_bps_h_0s = zip(concat_N_bps, concat_h_0s)
    concat_N_bps, concat_h_0s = zip(*sorted(zipped_N_bps_h_0s))
    concat_N_bps, concat_h_0s = list(concat_N_bps), list(concat_h_0s)

    #Curve fit split band/highest intensity band data with a degree 5 polynomial ( https://stackoverflow.com/questions/19165259/python-numpy-scipy-curve-fitting ):
    coeffs = np.polyfit(concat_N_bps, concat_h_0s, 5)
    poly_f = np.poly1d(coeffs)

    E_tots: List[float] = calc_E_tots(N_bps, Lks, poly_f, Lp, twist_Lp)
    unique_Lks: List[int] = sorted(set(Lks))
    E_tots_per_Lk: dict = calc_E_tots_per_Lk(N_bp_lim, unique_Lks, poly_f, Lp, twist_Lp)

    for N_bp, Lk in zip(N_bps, Lks):
        if Lk not in solid_range_per_Lk:
            solid_range_per_Lk[Lk] = [N_bp, N_bp]
        else:
            solid_range_per_Lk[Lk][0] = min(solid_range_per_Lk[Lk][0], N_bp)
            solid_range_per_Lk[Lk][1] = max(solid_range_per_Lk[Lk][1], N_bp)

    print("N_bp, E_tot, and calculated h_0 from polynomial fit (used to compute E_Tw) for lightgray curve:")
    for N_bp, E_tot in zip(N_bps, E_tots):
        print(N_bp, E_tot, poly_f(N_bp))
    print()

if not no_theory:
    N_bps_theory: List[int]; hs_theory: List[float]; Lks_theory: List[int]
    N_bps_theory, hs_theory, Lks_theory = calc_main_data(N_bp_lim_theory, turn_points_theory, Lk_lim_theory[0])

    E_tots_theory = calc_E_tots(N_bps_theory, Lks_theory, lambda x: h_canonical, Lp, twist_Lp)
    unique_Lks_theory: List[int] = sorted(set(Lks_theory))
    E_tots_per_Lk_theory: dict = calc_E_tots_per_Lk(N_bp_lim_theory, unique_Lks_theory, lambda x: h_canonical, Lp, twist_Lp)

    for N_bp, Lk in zip(N_bps_theory, Lks_theory):
        if Lk not in solid_range_per_Lk_theory:
            solid_range_per_Lk_theory[Lk] = [N_bp, N_bp]
        else:
            solid_range_per_Lk_theory[Lk][0] = min(solid_range_per_Lk_theory[Lk][0], N_bp)
            solid_range_per_Lk_theory[Lk][1] = max(solid_range_per_Lk_theory[Lk][1], N_bp)

    print("N_bp and total energy for fixed h_0 = 10.5 bp/turn (red data):")
    for N_bp, E_tot in zip(N_bps_theory, E_tots_theory):
        print(N_bp, E_tot)
    print()


apply_style()

#width_mm = 53
#height_mm = 37
#fig_width = width_mm / 25.4                 # in inches
#fig_height = height_mm / 25.4               # in inches

fig, ax = plt.subplots(figsize=fig_size_inches())

# major and minor tick locations
for axis in [ax.xaxis, ax.yaxis]:
    axis.set_minor_formatter(NullFormatter())

ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))
ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(10))

ax.yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(10))
ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))

ax.grid(which='major', axis='y', color='lightgray', linestyle='--', linewidth = LW.grid, zorder = -20, alpha=0.7) 
ax.grid(which='major', axis='x', color='lightgray', linestyle='--', linewidth = LW.grid, zorder = -20, alpha=0.7)

ax.set_xlim(*xlim)
ax.set_ylim(*ylim)

ax.set_xlabel(r'$N_{\mathregular{bp}}$', labelpad = 2)
ax.set_ylabel(r'$\Delta F/k_{\mathregular{B}}T$', labelpad = 2, rotation = 90)

#ax.set_position([0.125, 0.10999999999999999, 0.9, 0.88])

# plot
# dashed parabola extensions for each Lk, split into left and right segments to avoid a bridging line across the solid region (drawn first so solid lines appear on top)
if not only_theory:
    for Lk, (N_bps_full, E_tots_full) in E_tots_per_Lk.items():
        solid_min, solid_max = solid_range_per_Lk.get(Lk, [None, None])
        N_left  = [N for N in N_bps_full if N <= solid_min]
        E_left  = [E for N, E in zip(N_bps_full, E_tots_full) if N <= solid_min]
        N_right = [N for N in N_bps_full if N >= solid_max]
        E_right = [E for N, E in zip(N_bps_full, E_tots_full) if N >= solid_max]
        if N_left:
            ax.plot(N_left, E_left, '-', color='gray', linewidth=0.75, zorder=3)
            ax.plot(N_left, E_left, 'o', color='lightgray', markersize=MS.medium, markeredgecolor="none", zorder=4)
        if N_right:
            ax.plot(N_right, E_right, '-', color='gray', linewidth=0.75, zorder=3)
            ax.plot(N_right, E_right, 'o', color='lightgray', markersize=MS.medium, markeredgecolor="none", zorder=4)
    ax.plot(N_bps, E_tots, '-', color = 'gray', linewidth = 0.75, zorder = 5)          # lightgray solid lines
    ax.plot(N_bps, E_tots, 'o', color = 'lightgray', markersize = MS.medium, markeredgecolor="none", zorder = 5)        # lightgray solid markers
    special_x:      List[int]   = [N for N, Lk in zip(N_bps, Lks) if N in special_points.get(Lk, [])]
    special_E_tots: List[float] = [E for N, Lk, E in zip(N_bps, Lks, E_tots) if N in special_points.get(Lk, [])]
    hollow_x:       List[int]   = [N for N, Lk in zip(N_bps, Lks) if N in digested_points.get(Lk, [])]
    hollow_E_tots:  List[float] = [E for N, Lk, E in zip(N_bps, Lks, E_tots) if N in digested_points.get(Lk, [])]
    ax.plot(hollow_x, hollow_E_tots, 'o', color='white', markeredgecolor='black', markeredgewidth=LW.edge, markersize=MS.hollow, zorder=7)
    ax.plot(special_x, special_E_tots, 'o', color = 'black', markersize = MS.large, markeredgecolor="none", zorder = 6)

if not no_theory:
    for Lk, (N_bps_full, E_tots_full) in E_tots_per_Lk_theory.items():
        solid_min, solid_max = solid_range_per_Lk_theory.get(Lk, [None, None])
        N_left  = [N for N in N_bps_full if N <= solid_min]
        E_left  = [E for N, E in zip(N_bps_full, E_tots_full) if N <= solid_min]
        N_right = [N for N in N_bps_full if N >= solid_max]
        E_right = [E for N, E in zip(N_bps_full, E_tots_full) if N >= solid_max]
        if N_left:
            ax.plot(N_left, E_left, '-', color='gray', linewidth=0.75, zorder=3)
            ax.plot(N_left, E_left, 'o', color='lightgray', markersize=MS.medium, markeredgecolor="none", zorder=4)
        if N_right:
            ax.plot(N_right, E_right, '-', color='gray', linewidth=0.75, zorder=3)
            ax.plot(N_right, E_right, 'o', color='lightgray', markersize=MS.medium, markeredgecolor="none", zorder=4)
    ax.plot(N_bps_theory, E_tots_theory, '-', color = 'gray', linewidth = 0.75, zorder = 4)
    ax.plot(N_bps_theory, E_tots_theory, 'o', color = 'lightgray', markersize = MS.medium, markeredgecolor="none", zorder = 4)
    special_x:      List[int]   = [N for N, Lk in zip(N_bps_theory, Lks_theory) if N in special_points.get(Lk, [])]
    special_E_tots: List[float] = [E for N, Lk, E in zip(N_bps_theory, Lks_theory, E_tots_theory) if N in special_points.get(Lk, [])]
    hollow_x_theory:      List[int]   = [N for N, Lk in zip(N_bps_theory, Lks_theory) if N in digested_points_theory.get(Lk, [])]
    hollow_E_tots_theory: List[float] = [E for N, Lk, E in zip(N_bps_theory, Lks_theory, E_tots_theory) if N in digested_points_theory.get(Lk, [])]
    ax.plot(hollow_x_theory, hollow_E_tots_theory, 'o', color='white', markeredgecolor='black', markeredgewidth=LW.edge, markersize=MS.hollow, zorder=7)
    ax.plot(special_x, special_E_tots, 'o', color = 'black', markersize = MS.large, markeredgecolor="none", zorder = 5)

# tight layout to resize plot to fit plot + labels within canvas (bbox_inches = tight flag in savefig makes boundaries "stretch" to fit in the labels; this preserves cancas size as specified by fig_width and fig_height)
plt.tight_layout(pad = 0.2)

plot_base: str = plot_dir.rsplit('.', 1)[0] if '.' in plot_dir else plot_dir
plt.savefig(plot_base + '.png', dpi = 600)
plt.savefig(plot_base + '.svg', dpi = 600)

print(f"Saved plot to {plot_base}.png and {plot_base}.svg")
