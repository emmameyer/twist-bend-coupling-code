import sys
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator
from functions_for_plots import read_xvg_file
from fig_style import apply_style

"""
This script generates a stacking interaction colormap for DNA structures,
indicating intact and broken stacking interactions over time.
Initially written by Rachel Bricker, modified by Emma Meyer.

Plots base pairs on x-axis instead of individual residues.

A base pair is considered broken if either of the two paired residues has a broken stacking interaction.

example input:

python3 plot_stacking.py \
title \
com_files \
vec_files \
number of residues \
figure width \
figure height \
colorbar horizontal (1 for yes, 0 for vertical) \
topology (circular or linear) \
output directory \
path to scenarios

python3 plot_stacking.py \
"" \
com_files \
vec_files \
168 \
1.57480315 \
1.0 \
0 \
circular \
/Users/emma/Documents/kent/minicircles/OL21/analysis/D_84bp_deltaTw_minus_1_Tw_7_50ns_md/D_run1_84bp_deltaTw_minus_1_Tw_7_100ns_md/50ns/stacking_hbond \
/Users/emma/Documents/kent/minicircles/OL21/analysis/D_84bp_deltaTw_minus_1_Tw_7_50ns_md/D_run1_84bp_deltaTw_minus_1_Tw_7_100ns_md/50ns/stacking_hbond
"""

# system args
input_list      = sys.argv[1].replace("\\n", 
    "\n").replace("\\t", "\t").replace("\\(", "(").replace("\\)", ")")
legend          = input_list.split(',')
com_dir         = str(sys.argv[2])
vec_dir         = str(sys.argv[3])
n_residues      = int(sys.argv[4])
fig_width       = float(sys.argv[5])
fig_height      = float(sys.argv[6])
cbar_horizontal = bool(int(sys.argv[7]))
topology        = str(sys.argv[8]).lower()
output_dir_arg  = str(sys.argv[9])
paths           = list(sys.argv[10:])


# path and directory
for i in range(len(paths)):
    paths[i] = os.path.abspath(paths[i]) + "/"

output_dir = os.path.abspath(output_dir_arg) + "/"


# functions (defining the plane and such)
def S(alpha):
    return ( np.exp( -(alpha**4) ) +
             np.exp( -((alpha - np.pi)**4) ) +
             (0.1*np.exp( -((alpha - (0.5*np.pi))**4) )) )

def xi(COM_dist, alpha):
    return ( COM_dist / ( S(alpha) ) )


# processing functions
def get_data(paths, n_residues, com_dir, vec_dir, topology):
    stacking_coords = [ [] for _ in range(len(paths)) ]
    
    for scenario in range(len(paths)):
        time = []
        COM_coords = [ [] for _ in range(n_residues) ]
        norms = [ [] for _ in range(n_residues) ]

        # load COM coordinates
        for resi in range(n_residues):
            file = os.path.join(paths[scenario], com_dir,
                    f"nucleobase_COM_coord_{resi+1}.xvg")
            with open(file, "r") as f:
                data = read_xvg_file(f)
                for t in range(len(data)):
                    time_val = data[t].pop(0)
                    if resi == 0: time.append(time_val)
                    COM_coords[resi].append(data[t])

        # load VEC coords
        for resi in range(n_residues):
            file = os.path.join(paths[scenario], vec_dir, 
                    f"nucleobase_vec_coord_{resi+1}.xvg")
            with open(file, "r") as f:
                data = read_xvg_file(f)
                for t in range(len(data)):
                    data[t].pop(0)
                    atoms_xyz = [data[t][i:i+3] for i in range(0,
                                            len(data[t]), 3)]
                    COM_coord = COM_coords[resi][t]
                    vec_a = np.array(atoms_xyz[0]) - np.array(COM_coord)
                    vec_b = np.array(atoms_xyz[1]) - np.array(COM_coord)
                    norms[resi].append(np.cross(vec_a, vec_b))

        # stacking xi calculation
        for t in range(len(time)):
            time_step_xi = []
            mid = n_residues // 2

            # strand 1
            for i in range(0, mid):
                if topology == "circular":
                    next_i = (i + 1) % mid
                else:
                    if i == mid - 1: continue
                    next_i = i + 1

                p1, p2 = np.array(COM_coords[i][t]), np.array(COM_coords[next_i][t])
                dist = np.linalg.norm(p1 - p2)
                u, v = norms[i][t], norms[next_i][t]
                cos_theta = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
                alpha = np.arccos(np.clip(cos_theta, -1.0, 1.0))

                time_step_xi.append(xi(dist, alpha))
            
            # strand 2
            for i in range(mid, n_residues):
                if topology == "circular":                          # i added the circular logic
                    if i == n_residues - 1:
                        next_i = mid
                    else:
                        next_i = i + 1
                else:
                    if i == n_residues - 1: continue
                    next_i = i + 1
                
                p1, p2 = np.array(COM_coords[i][t]), np.array(COM_coords[next_i][t])
                dist = np.linalg.norm(p1 - p2)
                u, v = norms[i][t], norms[next_i][t]
                cos_theta = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
                alpha = np.arccos(np.clip(cos_theta, -1.0, 1.0))
                time_step_xi.append(xi(dist, alpha))
            
            stacking_coords[scenario].append(time_step_xi)

    return [t/1000 for t in time], stacking_coords


# convert residue to base pair
def convert_to_basepair_matrix(stacking_bool_matrix, n_residues):
    mid = n_residues // 2
    n_timepoints = stacking_bool_matrix.shape[0]
    
    basepair_matrix = np.zeros((n_timepoints, mid))
    
    for bp in range(mid):
        res_strand1 = bp
        res_strand2 = mid + (mid - bp - 1)
        
        stack_idx_strand1 = res_strand1
        stack_idx_strand2 = mid - 1 + (res_strand2 - mid + 1)
        
        for t in range(n_timepoints):
            broken_strand1 = stacking_bool_matrix[t, stack_idx_strand1] if stack_idx_strand1 < stacking_bool_matrix.shape[1] else 0.0
            broken_strand2 = stacking_bool_matrix[t, stack_idx_strand2] if stack_idx_strand2 < stacking_bool_matrix.shape[1] else 0.0
            
            basepair_matrix[t, bp] = max(broken_strand1, broken_strand2)
    
    return basepair_matrix


# plot
def plot_color_map_basepairs(time, stacking_bool_matrix, output_path):

    basepair_matrices = [convert_to_basepair_matrix(Z, n_residues) for Z in stacking_bool_matrix]

    fig, axes = plt.subplots(nrows=len(basepair_matrices), ncols=1,
                             sharex=True, figsize=(fig_width, fig_height))
    axes = np.atleast_1d(axes)
    cmap = mpl.colors.ListedColormap([(0.922, 0.922, 0.922), (0.62, 0.192, 0.961)])

    plt.subplots_adjust(hspace=0.025)

    time_arr = np.array(time)

    # Build time edges for pcolormesh: N timepoints -> N+1 edges
    dt = time_arr[1] - time_arr[0]
    time_edges = np.concatenate([[time_arr[0] - dt/2],
                                  (time_arr[:-1] + time_arr[1:]) / 2,
                                  [time_arr[-1] + dt/2]])

    # nick site separate plot, commented out for paper
    # nick_site_data = []
    # for Z_bp in basepair_matrices:
    #     nick_bp1 = Z_bp[:, 0:1]
    #     nick_bp2 = Z_bp[:, (mid_res-1):mid_res]
    #     nick_columns = np.concatenate([nick_bp1, nick_bp2], axis=1)
    #     nick_site_broken = np.any(nick_columns == 1, axis=1).astype(float)
    #     nick_site_data.append(nick_site_broken)

    # Nick site panel (top row) — disabled
    # nick_mesh = axes[0].pcolormesh(time_edges, [0, 1],
    #                                nick_site_data[0][np.newaxis, :],
    #                                cmap=cmap, vmin=0, vmax=1)
    # axes[0].set_ylabel("Nick\nsite", fontsize=font_size)
    # axes[0].set_yticks([0.5])
    # axes[0].set_yticklabels([''], fontsize=font_size-2)

    # Plot each scenario's base pair matrix
    for i, Z_bp in enumerate(basepair_matrices):
        n_bp = Z_bp.shape[1]

        bp_edges = np.arange(1, n_bp + 2)

        pixel_plot = axes[i].pcolormesh(time_edges, bp_edges,
                                        Z_bp.T,
                                        cmap=cmap, vmin=0, vmax=1)
        axes[i].set_title(legend[i])

        # Y-axis tick configuration (base pairs)
        y_major_freq = 20 if n_bp > 100 else 10
        y_minor_freq = 10 if n_bp > 100 else 5
        axes[i].set_yticks(np.arange(1, n_bp + 1, y_major_freq))
        axes[i].yaxis.set_minor_locator(FixedLocator(np.arange(1, n_bp + 1, y_minor_freq)))
        axes[i].tick_params(axis='y', which='minor', length=2, color='black')
        axes[i].tick_params(axis='both', which='major')

        x_minor_freq = 5
        x_major_freq = 10
        axes[i].set_xticks(np.arange(0, time_arr[-1] + x_major_freq, x_major_freq))
        axes[i].xaxis.set_minor_locator(FixedLocator(np.arange(0, time_arr[-1] + x_minor_freq, x_minor_freq)))
        axes[i].tick_params(axis='x', which='minor', length=2, color='black')
        axes[i].tick_params(axis='x', which='major')

        if i == len(basepair_matrices) // 2:
            axes[i].set_ylabel("Base Pair")

    axes[-1].set_xlabel("Time (ns)")

    if cbar_horizontal:
        cbar = fig.colorbar(pixel_plot, ax=axes, ticks=[], pad=0.3, orientation='horizontal')
        cbar.ax.text(0.35, 1.5, 'Intact', ha='center', va='center',
                      transform=cbar.ax.transAxes)
        cbar.ax.text(0.65, 1.5, 'Broken', ha='center', va='center',
                      transform=cbar.ax.transAxes)
    else:
        cbar = fig.colorbar(pixel_plot, ax=axes, ticks=[], pad=0.02, orientation='vertical')
        cbar.ax.text(2.0, 0.30, 'Intact', ha='center', va='center',
                     rotation=270, transform=cbar.ax.transAxes)
        cbar.ax.text(2.0, 0.70, 'Broken', ha='center', va='center',
                     rotation=270, transform=cbar.ax.transAxes)

    plt.savefig(output_path, bbox_inches='tight', dpi=600, pad_inches=0.02, facecolor='white')
    print(f"Base pair plot saved to: {output_path}")


# main
def main():
    time, coords = get_data(paths, n_residues, com_dir, vec_dir, topology)
    
    bool_mats = [(np.array(s) > 0.6).astype(float) for s in coords]
    apply_style()
    
    filename = f"stacking_map_basepairs_time_xaxis_{topology}.svg"
    full_output_path = os.path.join(output_dir, filename)
    
    plot_color_map_basepairs(time, bool_mats, full_output_path)

if __name__ == "__main__":
    main()
