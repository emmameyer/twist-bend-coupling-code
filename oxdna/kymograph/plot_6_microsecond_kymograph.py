"""
    Author: Rachel Bricker
    Year:   2023
"""

"""
    Plots a kymograph of the first 6 microseconds of a simulation.

    It is assumed there are enough configurations to plot 6 microseconds. msttcorefonts must be installed.
"""

"""
   usage: python3 plot_6_microsecond_kymograph.py
      1. [oxDNA input file]
      2. [path file to precise_kink_finder.py times output file]
      3. [path file to get_bend_and_twist_angles.py average bend and twist angle output file]
      4. [desired name of png file of kymograph (e.g., "six_microsecond_kymograph.png")]
"""

"""
    personal usage: 
        python3 /mnt/c/Users/brick/Documents/Programs/plot_kymograph/plot_6_microsecond_kymograph.py /mnt/c/Users/brick/Documents/Raw_data/84_circ_nick/input_84_nick /mnt/c/Users/brick/Documents/Raw_data/84_circ_nick/energetic_kink_finder_times_84_nick.txt /mnt/c/Users/brick/Documents/Raw_data/84_circ_nick/using_tommy_code/avg_twist_angle_data_2_prune.txt six_microsecond_kymograph_84_bp_ring.png
"""

import sys
from kymograph_functions import *

# command line inputs
input_file        = sys.argv[1]
stacked_times_txt = open(sys.argv[2], "r")
avg_angles_txt    = open(sys.argv[3], "r")
file_name         = str(sys.argv[4])

################################################################################################
#
# FUNCTIONS
#
################################################################################################

def plot_kymograph(colors : list[tuple[int]], avg_angles : list[float], file_name : str):
    """
        Creates a kymograph for the first 6 microseconds of a simulation

        Parameters:
            colors      (list[tuple[int]]) : list of rgb values (obtained from get_colors())
            avg_angles  (list[float])      : list of average twist angles (obtained from get_avg_angles())
            file_name   (str)              : name of the png file (must contain file type at the end
                                             such as ".svg" or ".png")
    """

    # fetch the number of configurations from colors list and average twist angle list
    n_confs = get_n_confs(colors, avg_angles)

    # create time (horizontal axis) list
    print_conf_interval, dt = get_print_conf_interval_and_dt(input_file)
    print_time_interval     = compute_print_time_interval(print_conf_interval, dt)    # units: nanoseconds
    print_time_interval     = print_time_interval / 1000    # change units to microseconds
    time                    = get_time_list(print_time_interval, n_confs)    # unit: microseconds (oxDNA time)

    # for font troubleshooting:
    # $ sudo apt install msttcorefonts -qq
    # $ rm ~/.cache/matplotlib -rf
    rcParameters = {'font.size': 7,                # font size
                    'font.family': "sans-serif",   # specify the font
                    'font.sans-serif': "Arial",
                    'axes.linewidth': 0.5,         # set width of axes
                    'xtick.major.width': 0.5,      # set tick widths
                    'xtick.minor.width': 0.5,
                    'ytick.major.width': 0.5,
                    'ytick.minor.width': 0.5}
    plt.rcParams.update(rcParameters)

    # the number of rows in the graph
    n_rows = 6

    # assign number of columns and rows
    fig, axes = plt.subplots(nrows=n_rows, ncols=1, gridspec_kw={'height_ratios': [1]*n_rows, 'width_ratios': [7]})

    for i in range(n_rows):
        # assign limits of horizontal axis
        x_lower_lim = i
        x_upper_lim = (i+1)

        # find where to slice lists
        index_1, index_2 = find_where_to_slice_list(time, x_lower_lim, x_upper_lim)

        # assign current axis
        ax1 = axes[i]

        # plot bp/turn
        base_pairs_per_turn = get_base_pairs_per_turn_list(avg_angles[index_1:index_2])

        # plot line to connect dots
        ax1.plot(time[index_1:index_2], base_pairs_per_turn, linewidth=0.5, color='lightgray', zorder=1)

        # plot dots
        ax1.scatter(time[index_1:index_2], base_pairs_per_turn, linewidth=0.5, c=colors[index_1:index_2], s=0.5, zorder=2)

        # secondary axis
        ax2 = ax1.secondary_yaxis('right', functions=(lambda x: 360/x, lambda x: 360/x))

        if i == (n_rows // 2):
            # label right y-axis for the middle row
            ax1.set_ylabel('$h$ (bp/turn)')

            # label left y-axis for the middle row
            ax2.set_ylabel('$\phi$ (${}^{\circ}$)', rotation=270, va='bottom')

        if (i+1) == len(axes):
            # label x-axis for only the bottom row
            plt.xlabel('Time ($\mu$s; oxDNA time)')

        # Set the range of x-axis
        ax1.set_xlim(x_lower_lim, x_upper_lim)

        # Set the range of y-axis s.t. 10.55 bp/turns is in the middle for all graphs
        ax1.set_ylim(10.55 - 0.75, 10.55 + 0.75)

        ax1.minorticks_on()
        ax2.minorticks_on()

    #fig.text(0.5, 0.015, 'Time ($\mu$s; oxDNA time)', ha='center')
    #fig.text(0.0075, 0.5, 'bp/turn', va='center', rotation='vertical')
    #fig.text(0.975, 0.5, r'$\langle\phi\rangle$ (degrees)', va='center', rotation='vertical')

    # figure dimensions
    set_size(7, 7)

    fig.tight_layout()

    # save figure
    plt.savefig(file_name, dpi=600, bbox_inches="tight")

################################################################################################
#
# MAIN PROGRAM
#
################################################################################################

def main():
    colors     = get_colors(stacked_times_txt)
    avg_angles = get_avg_angles(avg_angles_txt)

    plot_kymograph(colors, avg_angles, file_name)

if __name__ == "__main__": 
    main()
