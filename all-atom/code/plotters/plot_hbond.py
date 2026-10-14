"""
    Author: Rachel Bricker (and Emma Meyer)
    Year:   2025 (2026)
"""

"""
    Plots hydrogen bond analysis.
"""

"""
   usage: python3 plot_hbond.py
      1. list of model names for legend (must be parallel w.r.t the arguments i to i+n)
      2. name of .xvg file outputted by GROMACS utility `distance`
      3. name of .xvg file outputted by GROMACS utility `angle`
      4. enter 0 if constant temperature, 1 if simulated annealing
      5. figure width
      6. figure height
      7. horizontal color bar (1 for yes, 0 for no)
      8. output directory path
      9+. paths to data directories
   
   examples:
python3 plot_hbond_paper.py \
"" \
hbond_dist.xvg \
hbond_angle.xvg \
0 \
1.57480315 \
1.0 \
0 \
/Users/emma/Documents/kent/minicircles/OL21/analysis/D_84bp_deltaTw_minus_1_Tw_7_50ns_md/D_run1_84bp_deltaTw_minus_1_Tw_7_100ns_md/50ns/stacking_hbond \
/Users/emma/Documents/kent/minicircles/OL21/analysis/D_84bp_deltaTw_minus_1_Tw_7_50ns_md/D_run1_84bp_deltaTw_minus_1_Tw_7_100ns_md/50ns/stacking_hbond
"""

import sys
import numpy as np 
import statistics
from functions_for_plots import read_xvg_file
import matplotlib.pyplot as plt
from fig_style import apply_style
import matplotlib as mpl
from matplotlib.ticker import FixedLocator

# command line input
input_list      = sys.argv[1].replace("\\n", "\n").replace("\\t", "\t").replace("\\(", "(").replace("\\)", ")")
legend          = input_list.split(',')
dist_xvg        = str(sys.argv[2])
ang_xvg         = str(sys.argv[3])
annealing       = bool(int(sys.argv[4]))
fig_width       = float(sys.argv[5])
fig_height      = float(sys.argv[6])
cbar_horizontal = bool(int(sys.argv[7]))
output_dir_arg  = str(sys.argv[8])                                      # output directory path
paths           = list(sys.argv[9:])                                    # paths to data directories

# add trailing forward slash to directory path if necessary
for path in range(len(paths)):
    path_split = paths[path].split("/")
    if path_split[-1] != "":
        paths[path] = paths[path] + "/"

# set output directory
output_dir = output_dir_arg
if output_dir.split("/")[-1] != "":
    output_dir = output_dir + "/"

# functions

def get_dist(paths, dist_xvg):
    time      = []
    distances = [ [] for scenario in range(len(paths)) ]
    
    # loop over scenarios
    for scenario in range(len(paths)):
        file = paths[scenario] + dist_xvg
        with open(file, "r") as f:
            data = read_xvg_file(f) # read data file

            # loop over configurations
            for conf in range(len(data)):
                time_step = data[conf].pop(0)
                
                if scenario == 0:
                    time.append(time_step) # get time (measured in ps)

                # record distances
                distances[scenario].append(data[conf])
    
    time = [t/1000 for t in time] # convert ps -> ns
    
    return time, distances

def get_angle(paths, ang_xvg):
    angles = [ [] for scenario in range(len(paths)) ]
    
    # loop over scenarios
    for scenario in range(len(paths)):
        file = paths[scenario] + ang_xvg
        with open(file, "r") as f:
            data = read_xvg_file(f) # read data file

            # loop over configurations
            for conf in range(len(data)):
                data[conf].pop(0)

                # record distances
                angles[scenario].append(data[conf])
    
    return angles

def get_hbond_existence(distances, angles):
    # hbond exists if:
    #     distance <= 0.35 nm
    #     angle    <= 30 degrees

    # list of matrices
    hbond_bool_matrix = []

    for scenario in range(len(distances)):
        # initialize data matrix
        Z = []

        # number of configurations
        n_configurations = len(distances[scenario])

        # get number of base pairs
        base_pairs   = list(range(1, len(distances[scenario][0])+1))
        n_base_pairs = len(base_pairs)
        del base_pairs

        # loop over configurations
        for conf in range(n_configurations):
            if not len(Z): # if data matrix is empty
                Z = np.zeros((n_configurations, n_base_pairs)) # then create matrix of zeros

            # loop over base pairs
            for base_pair in range(n_base_pairs):
                dist  = distances[scenario][conf][base_pair]
                angle = angles[scenario][conf][base_pair]

                if not (dist <= 0.35 and angle <= 30):
                    Z[conf][base_pair]  = 1
        
        # append matrix
        hbond_bool_matrix.append(Z)

    return hbond_bool_matrix

def get_n_broken_hbond(hbond_bool_matrix, stop_residue_id=None):
    # hbond exists if:
    #     distance <= 0.35 nm
    #     angle    <= 30 degrees
    
    n_broken_hbond = [ [] for scenario in range(len(hbond_bool_matrix)) ]
    
    # loop over scenarios
    for scenario in range(len(hbond_bool_matrix)):

        # loop over configurations
        for conf in range(len(hbond_bool_matrix[scenario])):
            broken_hbond_count = 0
            
            if stop_residue_id == None:
                stop_residue_id = len(hbond_bool_matrix[scenario][conf])

            # loop over base pairs
            for base_pair in range(stop_residue_id):
                if hbond_bool_matrix[scenario][conf][base_pair] == 1:
                    broken_hbond_count += 1

            n_broken_hbond[scenario].append(broken_hbond_count)

    return n_broken_hbond
    
def get_avg_dist_per_conf(distances, stop_residue_id=None):  
    dist_avg = [ [] for scenario in range(len(distances)) ]
    
    # loop over scenarios
    for scenario in range(len(distances)):

        # loop over configurations
        for conf in range(len(distances[scenario])):
            
            if stop_residue_id == None:
                stop_residue_id = len(distances[scenario][conf])

            dist_avg[scenario].append(statistics.mean(distances[scenario][conf][:stop_residue_id]))
    
    return dist_avg

def plot_color_map(time, hbond_bool_matrix, annealing):
    font_size   = 6
    font_family = "Arial"

    file_name = output_dir + "colorplot_hbond_swapped_axes_binary.svg"
    if annealing:
        time = [t+600 for t in time]
        file_name = output_dir + "colorplot_hbond_binary_annealing.png"

    fig, axes = plt.subplots(nrows=len(hbond_bool_matrix), ncols=1,
                             sharex=True, figsize=(fig_width, fig_height))
    axes = np.atleast_1d(axes)

    cmap = mpl.colors.ListedColormap([(0.922, 0.922, 0.922), (0.62, 0.192, 0.961)])

    plt.subplots_adjust(hspace=0.025)

    time_arr = np.array(time)
    dt = time_arr[1] - time_arr[0]
    time_edges = np.concatenate([[time_arr[0] - dt/2],
                                  (time_arr[:-1] + time_arr[1:]) / 2,
                                  [time_arr[-1] + dt/2]])

    # nick site separate plot if wanted
    # nick_site_data = []
    # for Z in hbond_bool_matrix:
    #     n_base_pairs = Z.shape[1]
    #     nick_columns = np.concatenate([Z[:, 0:2], Z[:, n_base_pairs-2:n_base_pairs]], axis=1)
    #     nick_site_broken = np.any(nick_columns == 1, axis=1).astype(int)
    #     nick_site_broken = nick_site_broken[:, np.newaxis]
    #     nick_site_data.append(nick_site_broken)

    for i, Z in enumerate(hbond_bool_matrix):
        n_bp = Z.shape[1]

        # bp edges for pcolormesh y-axis
        bp_edges = np.arange(1, n_bp + 2)

        # Z is (n_timepoints, n_bp) -> .T is (n_bp, n_timepoints) = (M, N)
        pixel_plot = axes[i].pcolormesh(time_edges, bp_edges,
                                        Z.T,
                                        cmap=cmap, vmin=0, vmax=1)

        axes[i].set_title(legend[i])

        # y-axis ticks (base pairs)
        y_major_freq = 20 if n_bp > 100 else 10
        y_minor_freq = 10  if n_bp > 100 else 5
        axes[i].set_yticks(np.arange(1, n_bp + 1, y_major_freq))
        axes[i].yaxis.set_minor_locator(FixedLocator(np.arange(1, n_bp + 1, y_minor_freq)))
        axes[i].tick_params(axis='y', which='minor', length=2, color='black')

        x_major_freq = 10
        x_minor_freq = 5
        axes[i].set_xticks(np.arange(0, time[-1] + x_major_freq, x_major_freq))
        axes[i].xaxis.set_minor_locator(FixedLocator(np.arange(0, time[-1] + x_minor_freq, x_minor_freq)))
        axes[i].tick_params(axis='x', which='minor', length=2, color='black')

        if i == len(hbond_bool_matrix) // 2:
            axes[i].set_ylabel("Base Pair")

    axes[-1].set_xlabel("Time (ns)")


    if cbar_horizontal:
        cbar = fig.colorbar(pixel_plot, ax=axes, ticks=[], pad=0.3, orientation='horizontal')
        cbar.ax.text(0.35, 1.5, 'Intact', ha='center', va='center',
                     fontsize=font_size, transform=cbar.ax.transAxes)
        cbar.ax.text(0.65, 1.5, 'Broken', ha='center', va='center',
                     fontsize=font_size, transform=cbar.ax.transAxes)
    else:
        cbar = fig.colorbar(pixel_plot, ax=axes, ticks=[], pad=0.02, orientation='vertical')
        cbar.ax.text(2.0, 0.30, 'Intact', ha='center', va='center',
                     fontsize=font_size, rotation=270, transform=cbar.ax.transAxes)
        cbar.ax.text(2.0, 0.70, 'Broken', ha='center', va='center',
                     fontsize=font_size, rotation=270, transform=cbar.ax.transAxes)

    cbar.ax.tick_params(labelsize=font_size)
    for l in cbar.ax.yaxis.get_ticklabels():
        l.set_family(font_family)

    plt.savefig(file_name, bbox_inches="tight", dpi=600, pad_inches=0.02, facecolor='white')

# main

def main():  
    # get data from .xvg files
    time, distances   = get_dist(paths, dist_xvg)
    angles            = get_angle(paths, ang_xvg)
    hbond_bool_matrix = get_hbond_existence(distances, angles)
    n_broken_hbond    = get_n_broken_hbond(hbond_bool_matrix)

    # set rcParams
    apply_style()

    # plot boolean color map
    plot_color_map(time, hbond_bool_matrix, annealing)

    """
    # plot other data as function of time
    stop_residue_id = 6   # first six residues
    x_label         = "Simulation time (ns)"
    fig_width       = 3.35
    fig_height      = 1.4 
    plot_data(time,
              n_broken_hbond,
              x_label,
              "Number of broken\nbase pairs (bp)",
              None,
              legend,
              "melted_hbond_vs_time.svg",
              fig_width,
              fig_height)
    plot_data(time,
              get_avg_dist_per_conf(distances, stop_residue_id),
              x_label,
              r"$\langle r \rangle$ (nm)",
              None,
              legend,
              "distance_vs_time_terminal.svg",
              fig_width,
              fig_height)
    plot_data(time,
              get_n_broken_hbond(hbond_bool_matrix, stop_residue_id),
              x_label,
              "Number of broken base pairs (bp)",
              None,
              legend,
              "melted_hbond_vs_time_terminal.svg",
              fig_width,
              fig_height)
    """

    # print statistics
    frame_200ns  = 4000
    frame_1000ns = 20000
    for i in range(len(n_broken_hbond)):
        print("Average number of melted base pairs for file " + str(i+1) + ": " + str(round(statistics.mean(n_broken_hbond[i]),1)) + " +/- " + str(round(statistics.stdev(n_broken_hbond[i]),1)))
    for i in range(len(n_broken_hbond)):
        if len(n_broken_hbond[i]) > frame_200ns+1:
            print("Average number of melted base pairs for file " + str(i+1) + " (excluding first 200 ns): " + str(round(statistics.mean(n_broken_hbond[i][frame_200ns:]),1)) + " +/- " + str(round(statistics.stdev(n_broken_hbond[i][frame_200ns:]),1)))
    for i in range(len(n_broken_hbond)):
        if len(n_broken_hbond[i]) > frame_1000ns+1:
            print("Average number of melted base pairs for file " + str(i+1) + " (only including last 200 ns): " + str(round(statistics.mean(n_broken_hbond[i][frame_1000ns:]),1)) + " +/- " + str(round(statistics.stdev(n_broken_hbond[i][frame_1000ns:]),1)))

if __name__ == "__main__": 
    main()
