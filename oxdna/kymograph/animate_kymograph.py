"""
    Author: Rachel Bricker
    Year:   2023
"""

"""
   usage: python3 animate_kymograph.py
      1. [oxDNA input file]
      2. [path file to precise_kink_finder.py times output file]
      3. [path file to get_bend_and_twist_angles.py average bend and twist angle output file]
      4. [index of the configutation which you would like the movie to start at]
      5. [index of the configutation which you would like the movie to end at]
      6. [fps]
      7. [desired name of animation of kymograph (e.g., "kymograph_animation.mp4")]
"""

"""
    For personal usage:
    84 nick ring (typical print interval):
        python3 /mnt/c/Users/brick/Documents/Programs/plot_kymograph/animate_kymograph.py /mnt/c/Users/brick/Documents/Raw_data/84_circ_nick/input_84_nick /mnt/c/Users/brick/Documents/Raw_data/84_circ_nick/energetic_kink_finder_times_84_nick.txt /mnt/c/Users/brick/Documents/Raw_data/84_circ_nick/using_tommy_code/avg_twist_angle_data_2_prune.txt 10000 3610000 10 kymograph/kymograph_animation_84_nick.mp4
    84 nick ring (low print interval):
        python3 /mnt/c/Users/brick/Documents/Programs/plot_kymograph/animate_kymograph.py /mnt/c/Users/brick/Documents/Raw_data/84_circ_nick_low_print/input_84_nick_low_pint_align_trunc /mnt/c/Users/brick/Documents/Raw_data/84_circ_nick_low_print/energetic_kink_finder_times_84_nick_low_print.txt /mnt/c/Users/brick/Documents/Raw_data/84_circ_nick_low_print/using_tommy_code/avg_twist_angle_data_2_prune.txt 520800 2520800 24 kymograph/kymograph_animation_84_nick_low_pint.mp4
    89 nick ring (typical print interval):
        python3 /mnt/c/Users/brick/Documents/Programs/plot_kymograph/animate_kymograph.py /mnt/c/Users/brick/Documents/Raw_data/89_circ_nick/input_89_nick /mnt/c/Users/brick/Documents/Raw_data/89_circ_nick/energetic_kink_finder_times_89_nick.txt /mnt/c/Users/brick/Documents/Raw_data/89_circ_nick/using_tommy_code/avg_twist_angle_data_2_prune.txt 10000 3630000 10 kymograph/kymograph_animation_89_nick.mp4
    89 nick ring (low print interval)
        python3 /mnt/c/Users/brick/Documents/Programs/plot_kymograph/animate_kymograph.py /mnt/c/Users/brick/Documents/Raw_data/89_circ_nick_low_print/input_89_nick_low_pint /mnt/c/Users/brick/Documents/Raw_data/89_circ_nick_low_print/energetic_kink_finder_times_89_nick_low_pint.txt /mnt/c/Users/brick/Documents/Raw_data/89_circ_nick_low_print/using_tommy_code/avg_twist_angle_data_2_prune_low_pint.txt 400 345200 24 kymograph/kymograph_animation_89_nick_low_pint.mp4
"""

import sys
from matplotlib.collections import LineCollection
from matplotlib import animation
import numpy as np
from copy import deepcopy
from kymograph_functions import *

# command line inputs
input_file          = sys.argv[1]
stacked_times_txt   = open(sys.argv[2], "r")
avg_angles_txt      = open(sys.argv[3], "r")
starting_conf_index = int(sys.argv[4])
ending_conf_index   = int(sys.argv[5])
fps                 = int(sys.argv[6])
file_name           = str(sys.argv[7])

################################################################################################
#
# FUNCTIONS
#
################################################################################################

def make_rgb_paler(rgb : tuple[int]) -> tuple[int]:
    """
        Generates a rgb color code of a color which is a paler version of
        the given color. Function was obtained from
        https://stackoverflow.com/questions/33371939/calculate-rgb-equivalent-of-base-colors-with-alpha-of-0-5-over-white-background
        
        Parameters:
            rgb (tuple[int]) : the rgb color code of the color which
                               you would like to make paler
    """

    white = [1, 1, 1]
    alpha = 0.6
    return [alpha * c1 + (1 - alpha) * c2
            for (c1, c2) in zip(rgb, white)]

def get_moving_average(data_points, window_size=3) -> list[float]:
    """
        Calculates the moving (or rolling) mean of a set of data points. Reflective boundaries are
        assumed if the window size equals 3 such that the list returned has the same number of
        elements as the list passed. Function is based off of function from
        https://www.geeksforgeeks.org/how-to-calculate-moving-averages-in-python/

        Parameters:
            data_points     (list[int or float]) : data points to compute moving average from
            window_size     (int)                : the size of the window, i.e., the number of data
                                                   points used to calculate the moving average. 

        Returns:
            moving_averages (list[float])        : series of averages
    """

    moving_averages = []

    data_points_extended = deepcopy(data_points)

    if window_size == 3:
        # apply reflecting boundary conditions 
        data_points_extended.insert(0, data_points[0])
        data_points_extended.insert(-1, data_points[-1])

    for i in range(len(data_points_extended) - window_size + 1):
        # get the current window
        window = data_points_extended[i:(i + window_size)]

        # calculate the average of current window
        window_average = sum(window) / window_size

        # store the average of current window in moving average list
        moving_averages.append(window_average)

    return moving_averages

def init():
    """
        Function called once before the first frame to draw a clear frame.
    """

    solid_line.set_segments([])
    line_red.set_data([], [])
    line_blue.set_data([], [])

    return solid_line, line_red, line_blue,

def animate_kymograph(i):
    """
        Function to call at each frame.
    """

    # the below three lines of code is based of code provided by gg349 on Stack Overflow:
    # https://stackoverflow.com/questions/21077477/animate-a-line-with-different-colors
    # these three lines updates the mutlicolored solid line which follows the marker
    points   = np.array([time[:(i+1)], base_pair_per_turn_moving_avg[:(i+1)]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    solid_line.set_segments(segments)

    # the following if-else statement updates the position and color of the marker
    if paler_colors[i] == make_rgb_paler((1, 0, 0)):    # red marker
        line_red.set_data([time[i]], [base_pair_per_turn_moving_avg[i]])
        line_blue.set_data([], [])
    else:    # blue marker
        line_blue.set_data([time[i]], [base_pair_per_turn_moving_avg[i]])
        line_red.set_data([], [])

    return solid_line, line_red, line_blue,

################################################################################################
#
# MAIN PROGRAM
#
################################################################################################

print_conf_interval, dt = get_print_conf_interval_and_dt(input_file)
print_time_interval     = compute_print_time_interval(print_conf_interval, dt)    # units: nanoseconds
lower_bound             = convert_conf_index_to_time(starting_conf_index, print_time_interval, print_conf_interval)
upper_bound             = convert_conf_index_to_time(ending_conf_index, print_time_interval, print_conf_interval)
colors                  = get_colors(stacked_times_txt)
avg_angles              = get_avg_angles(avg_angles_txt)

# fetch the number of configurations from colors list and average twist angle list
n_confs = get_n_confs(colors, avg_angles)

# get times list (measured in nanoseconds hence print_time_interval is measured in nanoseconds)
time = get_time_list(print_time_interval, n_confs)

# shift time such that the graph starts at zero nanoseconds
time = [(t - lower_bound) for t in time]

# assign limits of horizontal axis such that the graph starts at zero nanoseonds
x_lower_lim = 0
x_upper_lim = (upper_bound - lower_bound)

# find where to slice lists to minimize size of svg file
index_1, index_2 = find_where_to_slice_list(time, x_lower_lim, x_upper_lim)

# slice lists
colors                        = colors[index_1:index_2]
avg_angles                    = avg_angles[index_1:index_2]
base_pair_per_turn            = get_base_pairs_per_turn_list(avg_angles)
base_pair_per_turn_moving_avg = get_moving_average(base_pair_per_turn)
time                          = time[index_1:index_2]

# set an appropriate figure width and marker size for the number of data points
if len(time) < 2500:
    fig_width   = 7
    marker_size = 7
else:
    fig_width   = 9
    marker_size = 8

# for font troubleshooting:
# $ sudo apt install msttcorefonts -qq
# $ rm ~/.cache/matplotlib -rf
rcParameters = {'font.size': 11,                # font size
                'font.family': "sans-serif",   # specify the font
                'font.sans-serif': "Arial",
                'axes.linewidth': 1,         # set width of axes
                'xtick.major.width': 1,      # set tick widths
                'xtick.minor.width': 1,
                'ytick.major.width': 1,
                'ytick.minor.width': 1}
plt.rcParams.update(rcParameters)

golden_mean = (np.sqrt(5)-1.0)/2.0     # aesthetic ratio
fig_height  = fig_width*golden_mean    # height in inches

fig, ax = plt.subplots(figsize=(fig_width, fig_height))

# plot grey line which connects the dots
ax.plot(time, base_pair_per_turn, linewidth=1, color=(0.91, 0.91, 0.91), zorder=1)

# make the color of the dots paler such that the marker in the animation is more visible
paler_colors = []
for rgb in colors:
    paler_colors.append(make_rgb_paler(rgb))

# plot dots
ax.scatter(time, base_pair_per_turn, linewidth=1, c=paler_colors, s=1, zorder=2)

# secondary axis
ax2 = ax.secondary_yaxis('right', functions=(lambda x: 360/x, lambda x: 360/x))

# for multicolored line
# https://stackoverflow.com/questions/21077477/animate-a-line-with-different-colors
# https://matplotlib.org/stable/gallery/shapes_and_collections/line_collection.html
solid_line = LineCollection([], colors=colors, linewidths=2, alpha=0.25)
solid_line.set_array(time)
ax.add_collection(solid_line)

# for marker (need two hence it can be either blue or red)
line_red,  = ax.plot([], [], 'o', markerfacecolor=(1, 0, 0), markersize=marker_size, markeredgecolor='w', markeredgewidth=1)
line_blue, = ax.plot([], [], 'o', markerfacecolor=(0, 0, 1), markersize=marker_size, markeredgecolor='w', markeredgewidth=1)

# label right y-axis
ax.set_ylabel('$h$ (bp/turn)')

# label x-axis
plt.xlabel('Time (ns; oxDNA time)', loc='center')

# label left y-axis
ax2.set_ylabel('$\phi$ (${}^{\circ}$)', rotation=270, va='bottom')

# Set the range of x-axis
ax.set_xlim(x_lower_lim, x_upper_lim)

# animate
FFwriter = animation.FFMpegWriter(fps=fps)
anim     = animation.FuncAnimation(fig, animate_kymograph, init_func=init, frames=len(time), interval=1, blit=False)

# print some information
exp_t = len(time)/fps
print("Expected video duration (s):", exp_t)
print("Number of confs:", len(time))

anim.save(file_name, writer=FFwriter, dpi=300)
