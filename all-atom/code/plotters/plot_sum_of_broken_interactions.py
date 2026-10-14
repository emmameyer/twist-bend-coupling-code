# this script finds the total number of broken hydrogen bonds 
# and the total number of broken stacking interactions for a circular dsDNA trajectory, 
# and plots the two values as a function of time.

"""
to use:

python3 plot_sum_of_broken_interactions.py \
--data_dir /Users/emma/Documents/kent/minicircles/OL21/analysis/I_84bp_deltaTw_minus_1_Lk_7_50ns_md/I_run1_84bp_deltaTw_minus_1_Lk_7_50ns_md/stacking_hbond \
--n_residues 168 \
--dist_xvg hbond_dist.xvg \
--ang_xvg hbond_angle.xvg \
--com_dir    com_files \
--vec_dir vec_files \
--title "" \
--output_dir /Users/emma/Documents/kent/sum_of_broken_interactions/I


required args:
--data_dir           path to directory (conatains all data files)
--n_residues         total number of residues
--dist_xvg           name of GROMACs distance .xvg file
--ang_xvg            name of GROMACs angle .xvg file
--com_dir            name of directory holding nucleobase COM .xvg files
--vec_dir            name of directory holding nucleobase VEC .xvg files
--output_dir         directory where output will be saved

optional args:
--title              plot title / legend label, default: ""
--fig_width          figure width in inches, default: 4.8
--fig_height         figure height in inches, default: 2.4

    [1] Jafilan, S., Klein, L., Hyun, C., & Florian, J. (2012). Intramolecular base stacking of dinucleoside monophosphate anions in aqueous solution.
        The Journal of Physical Chemistry B, 116(11), 3613-3618.
    [2] Condon, D. E., Kennedy, S. D., Mort, B. C., Kierzek, R., Yildirim, I., & Turner, D. H. (2015). Stacking in RNA: NMR of four tetramers benchmark molecular dynamics.
        Journal of chemical theory and computation, 11(6), 2729-2742.
"""

import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from functions_for_plots import read_xvg_file
from fig_style import apply_style, LW, TICK

# constants & parameters
stacking_threshold = 0.6                                        # value of xi, when xi > 0.6 stacking is broken
hbond_max_dist = 0.35                                           # nm       
hbond_max_angle = 30                                            # degrees

def _S(alpha):
    return ( np.exp( - ( alpha ** 4 )) + np.exp( - (( alpha - np.pi ) **4 )) 
        + 0.1 * np.exp( - (( alpha - 0.5 * np.pi ) ** 4 )) )

def _xi(com_dist, alpha):
    return com_dist / _S(alpha)

# load data

def load_xvg_time_and_columns(filepath):
    """
    Returns (time_ns list, list of row lists) from an .xvg file.
    """
    with open(filepath, "r") as f:
        raw = read_xvg_file(f)
    time_ns = []
    columns = []
    for row in raw:
        t = row[0]
        time_ns.append( t / 1000.0)
        columns.append(row[1:])
    return time_ns, columns

def get_hbond_data(data_dir, dist_xvg, ang_xvg):
    """
    Returns: 
    time_ns: list of float; and n_broken: list of int, 
    the number of broken hbonds at each time step.
    """
    dist_path = os.path.join(data_dir, dist_xvg)
    ang_path = os.path.join(data_dir, ang_xvg)

    time_ns, distances = load_xvg_time_and_columns(dist_path)
    _,       angles    = load_xvg_time_and_columns(ang_path)

    n_broken = []
    for conf in range(len(distances)):
        count = 0
        for bp in range(len(distances[conf])):
            dist = distances[conf][bp]
            angle = angles[conf][bp]
            if not (dist <= hbond_max_dist and angle <= hbond_max_angle):
                count += 1
        n_broken.append(count)
    
    return time_ns, n_broken

def get_stacking_data(data_dir, n_residues, com_dir, vec_dir):
    """
    Returns:
    time_ns: list of float; n_broken: list of int
    number of base pairs with at least one broken stacking interaction at each time step.
    """
    mid = n_residues // 2

    com_path = os.path.join(data_dir, com_dir)
    vec_path = os.path.join(data_dir, vec_dir)

    # com coords
    COM_coords = []
    time_ns = None
    for resi in range(n_residues):
        fpath = os.path.join(com_path, f"nucleobase_COM_coord_{resi+1}.xvg")
        t, cols = load_xvg_time_and_columns(fpath)
        if time_ns is None:
            time_ns = t
        COM_coords.append(cols)

    # vec coords
    norms = []
    for resi in range(n_residues):
        fpath = os.path.join(vec_path, f"nucleobase_vec_coord_{resi+1}.xvg")
        _, cols = load_xvg_time_and_columns(fpath)
        resi_norms = []
        for t_idx, row in enumerate(cols):
            atoms_xyz = [row[i:i+3] for i in range (0, len(row), 3)]
            com       = np.array(COM_coords[resi][t_idx])
            vec_a     = np.array(atoms_xyz[0]) - com
            vec_b     = np.array(atoms_xyz[1]) - com
            resi_norms.append(np.cross(vec_a, vec_b))
        norms.append(resi_norms)

    # compute xi for each stacking pair
    n_frames = len(time_ns)
    stacking_xi = []
    for t in range(n_frames):
        frame_xi = []

        # strand 1
        for i in range(mid):
            next_i = (i+1) % mid
            p1 = np.array(COM_coords[i][t])
            p2 = np.array(COM_coords[next_i][t])
            dist = np.linalg.norm(p1 - p2)
            u, v = norms[i][t], norms[next_i][t]
            cos_theta = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
            alpha = np.arccos(np.clip(cos_theta, -1.0, 1.0))
            frame_xi.append(_xi(dist, alpha))

        # strand 2
        for i in range(mid, n_residues):
            next_i = mid if i == n_residues -1 else i +1
            p1 = np.array(COM_coords[i][t])
            p2 = np.array(COM_coords[next_i][t])
            dist = np.linalg.norm(p1 - p2)
            u, v = norms[i][t], norms[next_i][t]
            cos_theta = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
            alpha = np.arccos(np.clip(cos_theta, -1.0, 1.0))
            frame_xi.append(_xi(dist, alpha))

        stacking_xi.append(frame_xi)

    # convert to boolean matrix
    xi_array = np.array(stacking_xi)
    bool_mat = (xi_array > stacking_threshold).astype(float)

    # convert to base pairs. a bp is broken if either two residues has a broken
    # stacking interaction.
    n_broken = []
    for t in range(n_frames):
        count = 0 
        for bp in range(mid):
            res_strand1         = bp
            res_strand2         = mid + (mid - bp - 1)
            stack_idx_strand1   = res_strand1
            stack_idx_strand2   = mid -1 + (res_strand2 - mid + 1)

            broken_s1 = bool_mat[t, stack_idx_strand1] if stack_idx_strand1 < bool_mat.shape[1] else 0.0
            broken_s2 = bool_mat[t, stack_idx_strand2] if stack_idx_strand2 < bool_mat.shape[1] else 0.0

            if max(broken_s1, broken_s2) == 1.0:
                count += 1
        n_broken.append(count)

    return time_ns, n_broken

# rolling avg
def smooth(y, window=50):
    y = np.array(y)
    pad = window // 2
    y_padded = np.pad(y, (pad, pad), mode='edge')
    smoothed = np.convolve(y_padded, np.ones(window)/window, mode='valid')
    return smoothed[:len(y)]

# plotting
def plot_single_axis(time_ns, n_broken_hbond, n_broken_stack, title, output_path):
    apply_style()
    font_size = 6

    width_mm = 28
    height_mm = 24
    width = width_mm / 25.4
    height = height_mm / 25.4

    fig, ax = plt.subplots(figsize=(width, height))

    color_hbond = (1.0, 0.0, 0.0) 
    color_stack = (0.0, 0.0, 1.0)

    #ax.plot(time_ns, n_broken_hbond, color=color_hbond, linewidth=0.8, label="Broken H-bonds")    rough data
    #ax.plot(time_ns, n_broken_stack, color=color_stack, linewidth=0.8, label="Broken Stacking")

    ax.plot(time_ns, smooth(n_broken_hbond, window=50), color=color_hbond, linewidth=LW.data, label="H-bonds")
    ax.plot(time_ns, smooth(n_broken_stack, window=50), color=color_stack, linewidth=LW.data, label="Stacking")

    ax.set_xlim(0, time_ns[-1])
    ax.set_ylim(0, max(max(n_broken_hbond[1:]), max(n_broken_stack[1:])) + 1)
    ax.set_xlabel("", labelpad = TICK.pad_major)#Time (ns)", fontsize=font_size)
    ax.set_ylabel("", labelpad = TICK.pad_major)## Broken", fontsize=font_size)

    ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))
    ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(1))
    ax.tick_params(axis='both', labelsize = font_size, pad = TICK.pad_major)


    if title:
        ax.set_title(title, fontsize=font_size)
    
    ax.legend(fontsize=font_size -1, frameon=False, loc='upper right')

    #fig.tight_layout()
    plt.tight_layout(pad = 0.3)
    plt.savefig(output_path, bbox_inches="tight", dpi=300, facecolor='white')
    print(f"Saved plot to {output_path}")

    #print diagnostics
    #print average number of broken hbonds/stacking over the whole trajectory, excluding the first frame
    avg_hbonds = np.mean(n_broken_hbond[1:])
    avg_stacking = np.mean(n_broken_stack[1:])
    print(f"Average number of broken hbonds (excluding first frame): {avg_hbonds:.2f}")
    print(f"Average number of broken stacking interactions (excluding first frame): {avg_stacking:.2f}")


# parsing args
def parse_args():
    parser = argparse.ArgumentParser(description="Plot number of broken hbond and stacking interactions over time for a circular dsDNA trajectory.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    req = parser.add_argument_group("required arguments")
    req.add_argument("--data_dir", required=True, help="Path to directory containing all data files.")
    req.add_argument("--n_residues", required=True, type=int, help="Total number of residues.")
    req.add_argument("--dist_xvg", required=True, help="Name of gmx distance .xvg file.")
    req.add_argument("--ang_xvg", required=True, help="Name of gmx angle .xvg file.")
    req.add_argument("--com_dir", required=True, help="Name of directory holding COM .xvg files.")
    req.add_argument("--vec_dir", required=True, help="Name of directory holding VEC .xvg files.")
    req.add_argument("--output_dir", required=True, help="Directory where you want output.")
    
    parser.add_argument("--title", default="", help="Plot title / legend label.")
    parser.add_argument("--fig_width", default=4.8, type=float, help="Figure width in inches.")
    parser.add_argument("--fig_height", default=2.4, type=float, help="Figure height in inches.")

    return parser.parse_args()


# main
def main():
    data_dir = os.path.abspath(args.data_dir)
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)

    print("Loading hbond data")
    time_ns, n_broken_hbond = get_hbond_data(data_dir, args.dist_xvg, args.ang_xvg)

    print("Loading stacking data")
    time_ns_stack, n_broken_stack = get_stacking_data(data_dir, args.n_residues, args.com_dir, args.vec_dir)

    output_path = os.path.join(output_dir, "broken_hbond_stacking_over_time.svg")
    plot_single_axis(time_ns, n_broken_hbond, n_broken_stack, args.title, output_path)

if __name__ == "__main__":
    args = parse_args()
    main()
