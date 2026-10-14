"""
    Author: Rachel Bricker
    Year:   2023
"""

"""
    Plots a 500 ns kymograph. It is assumed there are enough configurations to plot 500 ns after
    the starting configuration provided. msttcorefonts must be installed.
"""

"""
   usage: python3 plot_500_ns_kymograph.py
      1. [path to oxDNA input file]
      2. [path to precise_kink_finder.py times output file]
      3. [path to get_bend_and_twist_angles.py average twist angle output file]
      4. [index of the configutation which you would like the graph to start at]
      5. [show left y-axis (enter 1 for yes, 0 for no)]
      6. [label the x-axis (enter 1 for yes, 0 for no)]
      7. [show secondary y-axis (enter 1 for yes, 0 for no)]
      8. [desired name of kymograph (e.g., "kymograph_82_nick_circ.svg")]
"""

"""
    For personal usage:
    82 bp ring:
        python3 /mnt/c/Users/brick/Documents/Programs/plot_kymograph/plot_500_ns_kymograph.py /mnt/c/Users/brick/Documents/Raw_data/82_circ_nick/input_82_nick /mnt/c/Users/brick/Documents/Raw_data/82_circ_nick/energetic_kink_finder_times_82_nick.txt /mnt/c/Users/brick/Documents/Raw_data/82_circ_nick/using_tommy_code/avg_twist_angle_data_2_prune.txt 106660000 1 0 0 kymograph_82_nick_circ.svg
    84 bp ring:
        python3 /mnt/c/Users/brick/Documents/Programs/plot_kymograph/plot_500_ns_kymograph.py /mnt/c/Users/brick/Documents/Raw_data/84_circ_nick/input_84_nick /mnt/c/Users/brick/Documents/Raw_data/84_circ_nick/energetic_kink_finder_times_84_nick.txt /mnt/c/Users/brick/Documents/Raw_data/84_circ_nick/using_tommy_code/avg_twist_angle_data_2_prune.txt 46660000 0 1 0 kymograph_84_nick_circ.svg
    86 bp ring:
        python3 /mnt/c/Users/brick/Documents/Programs/plot_kymograph/plot_500_ns_kymograph.py /mnt/c/Users/brick/Documents/Raw_data/86_circ_nick/input_86_nick /mnt/c/Users/brick/Documents/Raw_data/86_circ_nick/energetic_kink_finder_times_86_nick.txt /mnt/c/Users/brick/Documents/Raw_data/86_circ_nick/using_tommy_code/avg_twist_angle_data_2_prune.txt 80000000 0 0 1 kymograph_86_nick_circ.svg
"""

import sys
from kymograph_functions import *

# command line inputs
input_file            = sys.argv[1]
stacked_times_txt     = open(sys.argv[2], "r")
avg_angles_txt        = open(sys.argv[3], "r")
starting_conf_index   = int(sys.argv[4])
show_left_y_axis      = int(sys.argv[5])
label_x_axis          = int(sys.argv[6])
show_secondary_y_axis = int(sys.argv[7])
file_name             = str(sys.argv[8])

################################################################################################
#
# FUNCTIONS
#
################################################################################################

def plot_kymograph(colors : list[tuple[int]], avg_angles : list[float], starting_conf_index : int, file_name : str,
                   show_left_y_axis : int, label_x_axis : int, show_secondary_y_axis : int):
    """
        Creates 500 ns kymograph.

        Parameters:
            colors                (list[tuple[int]]) : list of rgb values (obtained from get_colors())
            avg_angles            (list[float])      : list of average twist angles (obtained from get_avg_angles())
            starting_conf_index   (int)              : index of configuration where graph begins
            file_name             (str)              : name of the png file (must contain file type at the end
                                                       such as ".svg" or ".png")
            show_left_y_axis      (int)              : flags whether the main (left) y-axis should be shown or not
            label_x_axis          (int)              : flags whether the x-axis should be labeled or not
            show_secondary_y_axis (int)              : flags whether the secondary (right) y-axis should be shown or not
    """

    # fetch the number of configurations from colors list and average twist angle list
    n_confs = get_n_confs(colors, avg_angles)

    # create time (horizontal axis) list
    print_conf_interval, dt = get_print_conf_interval_and_dt(input_file)
    print_time_interval     = compute_print_time_interval(print_conf_interval, dt)    # units: nanoseconds
    time                    = get_time_list(print_time_interval, n_confs)    # unit: microseconds (oxDNA time)

    # get lower bound in ns
    lower_bound = convert_conf_index_to_time(starting_conf_index, print_time_interval, print_conf_interval)
    
    # assign time upper bound such that the range is 500 nanoseconds
    upper_bound = lower_bound + 500

    # shift time such that the graph starts at zero nanoseconds
    time = [(t - lower_bound) for t in time]

    # assign limits of horizontal axis
    x_lower_lim = 0
    x_upper_lim = (upper_bound - lower_bound)

    # find where to slice lists to minimize size of svg file
    index_1, index_2 = find_where_to_slice_list(time, x_lower_lim, x_upper_lim)

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

    fig, ax1 = plt.subplots()

    # get bp/turn
    base_pairs_per_turn = get_base_pairs_per_turn_list(avg_angles[index_1:index_2])

    # plot line to connect dots
    ax1.plot(time[index_1:index_2], base_pairs_per_turn, linewidth=0.5, color='lightgray', zorder=1)

    # plot dots
    ax1.scatter(time[index_1:index_2], base_pairs_per_turn, linewidth=0.5, c=colors[index_1:index_2], s=0.5, zorder=2)

    # create secondary axis
    ax2 = ax1.secondary_yaxis('right', functions=(lambda x: 360/x, lambda x: 360/x))

    if show_left_y_axis:
        # label left y-axis
        ax1.set_ylabel('$h$ (bp/turn)')
    else:
        # remove axis values
        ax1.set_yticks([])

    if show_secondary_y_axis:
        # label right y-axis
        ax2.set_ylabel('$\phi$ (${}^{\circ}$)', rotation=270, va='bottom')
    else:
        # remove axis values
        ax2.set_yticks([])

    if label_x_axis:
        # label x-axis
        plt.xlabel('Time (ns; oxDNA time)')

    # set the range of x-axis
    ax1.set_xlim(x_lower_lim, x_upper_lim)

    # set the range of y-axis s.t. 10.55 bp/turns is in the middle for all graphs
    ax1.set_ylim(10.55 - 0.55, 10.55 + 0.55)

    # make minor ticks visible
    ax1.minorticks_on()
    ax2.minorticks_on()

    # figure dimensions
    set_size(2.3, 1)

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

    plot_kymograph(colors, avg_angles, starting_conf_index, file_name,
                   show_left_y_axis, label_x_axis, show_secondary_y_axis)

if __name__ == "__main__": 
    main()
