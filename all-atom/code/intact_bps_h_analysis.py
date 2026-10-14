# python script to analyze the helical repeat of the intact base pairs for the last 10 ns of a simulation

"""
intact_bps_h_analysis.py

finds the avergae helical repeat (bp/turn) of only the intact base pairs of a circular dsDNA, averaged over just the last 10 ns

base pairs are considered intact only if both of these conditions are true:
1. h bond is intact when: dist <= 0.35 nm and angle <= 30 degrees for donor-acceptor-hydrogen
2. stacking is intact when: xi <= 0.6 for both residues of the base pair

A twist step between bp i and bp i+1 is included only if both of them are intact

to use:

python3 intact_bps_h_analysis.py \
--data_dir /Users/emma/Documents/kent/projects/minicircles/OL21/analysis/I_84bp_deltaTw_minus_1_Lk_7_50ns_md/I_run1_84bp_deltaTw_minus_1_Lk_7_50ns_md/stacking_hbond \
--n_residues 168 \
--n_valid_steps 84 \
--dist_xvg hbond_dist.xvg \
--ang_xvg hbond_angle.xvg \
--com_dir com_files \
--vec_dir vec_files \
--twist_xvg /Users/emma/Documents/kent/projects/minicircles/OL21/analysis/I_84bp_deltaTw_minus_1_Lk_7_50ns_md/I_run1_84bp_deltaTw_minus_1_Lk_7_50ns_md/twist/Twist__twist.xvg \
--last_ns 10

required arguments:
--data_dir        directory that has stacking/hbond data
--n_residues      number of residues
--n_valid_steps   number of valid sequential twist steps in do_x3dna output (for circular dna it should be Nbp)
--dist_xvg        name of GROMACS distance .xvg file
--ang_xvg         name of GROMACS angle .xvg file
--com_dir         name of subdirectory holding nucleobase COM .xvg files
--vec_dir         name of subdirectory holding nucleobase vector .xvg files
--twist_xvg       path to Twist__twist.xvg from do_x3dna

optional arguments:
--last_ns         how many nanoseconds from the end of the simulation to analyze (default: 10)
--hbond_max_dist  h bond distance cutoff in nm (default: 0.35)
--hbond_max_angle h bond angle cutoff in degrees (default: 30)
--stacking_thresh xi threshold above which stacking is broken (default: 0.6)
"""

import os
import argparse
import numpy as np

# stacking parameters

def _S(alpha):
    """Orientational factor used in the xi stacking metric."""
    return (np.exp(-(alpha ** 4))
            + np.exp(-((alpha - np.pi) ** 4))
            + 0.1 * np.exp(-((alpha - 0.5 * np.pi) ** 4)))


def _xi(com_dist, alpha):
    return com_dist / _S(alpha)


# load files

def read_xvg(filepath):
    """
    Parse an .xvg file, skipping lines that start with # or @.
    Returns a 2-D numpy array: rows = frames, columns = all values including time.
    """
    rows = []
    with open(filepath, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line[0] in ("#", "@"):
                continue
            rows.append([float(v) for v in line.split()])
    return np.array(rows)


# determine where there are h bonds and if its intact

def get_hbond_intact(data_dir, dist_xvg, ang_xvg,
                     hbond_max_dist, hbond_max_angle):
    """
    Returns
    -------
    time_ps : np.ndarray, shape (n_frames,)
    intact  : np.ndarray of bool, shape (n_frames, n_bp)
        intact[t, bp] is True when the H-bond of base pair bp is intact at
        frame t.
    """
    dist_data = read_xvg(os.path.join(data_dir, dist_xvg))
    ang_data  = read_xvg(os.path.join(data_dir, ang_xvg))

    time_ps   = dist_data[:, 0]
    distances = dist_data[:, 1:]   # shape (n_frames, n_bp)
    # Column 1 of the angle file is the per-frame average written by GROMACS
    # gmx angle; skip it so the remaining columns align with the distance data.
    angles    = ang_data[:, 2:]    # shape (n_frames, n_bp)

    intact = (distances <= hbond_max_dist) & (angles <= hbond_max_angle)
    return time_ps, intact


# determine where there is stacking and if its intact

def get_stacking_intact(data_dir, n_residues, com_dir, vec_dir,
                        stacking_thresh):
    """
    Returns
    -------
    time_ps      : np.ndarray, shape (n_frames,)
    bp_intact    : np.ndarray of bool, shape (n_frames, n_bp)
        bp_intact[t, bp] is True when BOTH stacking interactions that involve
        base pair bp are intact at frame t (i.e. the stacking step on each
        side of that bp does not exceed the threshold).

    For a circular molecule with n_bp base pairs:
      - Strand 1 residues: 0 … n_bp-1  (0-indexed, so 0 here is bp 1, and bp-1 is n_bp)
      - Strand 2 residues: n_bp … 2*n_bp-1

    Stacking index on strand 1:
      step i connects residue i to residue (i+1) % n_bp: 84 steps total

    Stacking index on strand 2:
      step i (0-indexed within strand 2) connects residue (n_bp + i) to residue (n_bp + (i+1) % n_bp)

    Base-pair mapping (antiparallel):
      bp index j  ->  strand-1 residue j,  strand-2 residue (n_bp - 1 - j) + n_bp
      (residue j pairs with residue 2*n_bp - 1 - j)

    A bp is stacking-intact if all four stacking steps that touch its two
    residues are intact.  In practice we check the two steps on each strand
    that border the bp's residue:
      strand-1 step to the left of residue j:  (j - 1) % n_bp
      strand-1 step to the right of residue j: j % n_bp
    (and analogously for the strand-2 residue).

    --  --          / --        if the four residues (shown on the left, straight and fine) are intact, then that twist step between them contributes to the "intact" h
    --  --         -- /         if the four residues (shown on the right, flipped out and broken) are not intact, then that twist step does not contribute to the "intact" h
    """
    n_bp  = n_residues // 2
    n_res = n_residues

    com_path = os.path.join(data_dir, com_dir)
    vec_path = os.path.join(data_dir, vec_dir)

    # Load COM coordinates ─ shape per residue: (n_frames, 3)
    COM_coords = []
    time_ps    = None
    for resi in range(n_res):
        d = read_xvg(os.path.join(com_path, f"nucleobase_COM_coord_{resi+1}.xvg"))
        if time_ps is None:
            time_ps = d[:, 0]
        COM_coords.append(d[:, 1:4])   # x, y, z

    # Load VEC coordinates and compute normal vectors ─ shape per residue: (n_frames, 3)
    norms = []
    for resi in range(n_res):
        d = read_xvg(os.path.join(vec_path, f"nucleobase_vec_coord_{resi+1}.xvg"))
        cols = d[:, 1:]                # everything after time
        com  = COM_coords[resi]        # (n_frames, 3)
        # Each row has two atoms: [x1 y1 z1 x2 y2 z2]
        atom_a = cols[:, 0:3] - com
        atom_b = cols[:, 3:6] - com
        resi_norms = np.cross(atom_a, atom_b)  # (n_frames, 3)
        norms.append(resi_norms)

    n_frames = len(time_ps)

    # Compute xi for every stacking step on both strands
    # xi_strand1[t, i] = xi for step i on strand 1  (i: 0 … n_bp-1)
    # xi_strand2[t, i] = xi for step i on strand 2

    def compute_xi_strand(res_indices, wrap):
        """
        res_indices : list of 0-based residue indices for this strand, in order
        wrap        : if True, last step wraps back to first residue (circular)
        """
        n = len(res_indices)
        xi_mat = np.zeros((n_frames, n))
        for i in range(n):
            next_i = (i + 1) % n if wrap else (i + 1)
            if not wrap and next_i >= n:
                break
            r1 = res_indices[i]
            r2 = res_indices[next_i]
            p1 = COM_coords[r1]   # (n_frames, 3)
            p2 = COM_coords[r2]
            dist = np.linalg.norm(p1 - p2, axis=1)  # (n_frames,)
            u = norms[r1]
            v = norms[r2]
            cos_theta = np.einsum('ij,ij->i', u, v) / (
                np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1))
            alpha = np.arccos(np.clip(cos_theta, -1.0, 1.0))
            xi_mat[:, i] = _xi(dist, alpha)
        return xi_mat

    strand1_res = list(range(n_bp))
    strand2_res = list(range(n_bp, n_res))

    xi_s1 = compute_xi_strand(strand1_res, wrap=True)   # (n_frames, n_bp)
    xi_s2 = compute_xi_strand(strand2_res, wrap=True)   # (n_frames, n_bp)

    broken_s1 = xi_s1 > stacking_thresh   # (n_frames, n_bp)  step i on strand 1
    broken_s2 = xi_s2 > stacking_thresh

    # For each bp j, determine if its stacking is intact.
    # Strand-1 residue for bp j is index j.
    # Strand-2 residue for bp j (antiparallel) is index n_bp - 1 - j on strand 2
    # (i.e., absolute residue index n_bp + (n_bp - 1 - j)).
    # A stacking step on strand 1 is broken if broken_s1[:, step_i] is True.
    # Step i connects residue i to residue (i+1) % n_bp.
    # Residue j on strand 1 is touched by steps (j-1)%n_bp (left) and j (right).
    # Same logic for strand 2 with the reversed indexing.

    bp_intact_stack = np.ones((n_frames, n_bp), dtype=bool)
    for j in range(n_bp):
        step_left_s1  = (j - 1) % n_bp
        step_right_s1 = j

        s2_res_local   = n_bp - 1 - j          # local index within strand 2
        step_left_s2   = (s2_res_local - 1) % n_bp
        step_right_s2  = s2_res_local

        bp_intact_stack[:, j] = (
            ~broken_s1[:, step_left_s1]  &
            ~broken_s1[:, step_right_s1] &
            ~broken_s2[:, step_left_s2]  &
            ~broken_s2[:, step_right_s2]
        )

    return time_ps, bp_intact_stack


# load twist data from do_x3dna output

def load_twist(twist_xvg, n_valid_steps):
    """
    Load twist data from a do_x3dna .xvg file.

    parameters:
    twist_xvg     : path to the file
    n_valid_steps : number of genuine sequential base-pair steps, (columns 1 … n_valid_steps in the data; everything beyond that is a fantasized step by do_x3dna and is discarded)

    returns:
    time_ps : np.ndarray, shape (n_frames,)
    twist   : np.ndarray, shape (n_frames, n_valid_steps)
              twist[t, i] = twist of step i at frame t, in degrees.

    for a circular N-bp molecule the valid steps are 0-indexed 0 … N-1.
    step i connects bp i to bp (i+1) % N.
    """
    data = read_xvg(twist_xvg)
    time_ps = data[:, 0]
    twist   = data[:, 1 : n_valid_steps + 1]
    return time_ps, twist


# find h

def compute_intact_helical_repeat(time_ps_hb, hb_intact,
                                  time_ps_st, st_intact,
                                  time_ps_tw, twist,
                                  last_ns, n_valid_steps):
    """
    for each frame in the last 'last_ns' nanoseconds:
      1. determine which base pairs are fully intact
      2. identify twist steps where all four residues of the two bps are intact
      3. collect the twist values for those steps (filtering out invalid do_x3dna values)
      4. average over all valid intact steps in the frame

    output:
      - per-frame averages
      - grand mean and block-averaging standard error across all frames

    returns:
      dict with results
    """
    last_ps = last_ns * 1000.0

    # Align frames: all three datasets should share the same time axis.
    time_max    = time_ps_hb[-1]
    cutoff_ps   = time_max - last_ps

    twist_time_to_idx = {round(t, 4): i for i, t in enumerate(time_ps_tw)}

    frame_avgs   = []
    frame_counts = []
    included_frames = 0

    n_bp = hb_intact.shape[1]

    for t_idx, t in enumerate(time_ps_hb):
        if t <= cutoff_ps:
            continue

        # find the corresponding twist frame
        t_key = round(t, 4)
        if t_key not in twist_time_to_idx:
            continue

        tw_idx = twist_time_to_idx[t_key]

        # fully intact base pairs at this frame
        bp_intact = hb_intact[t_idx] & st_intact[t_idx]

        # collect twist values for steps where all four residues are intact
        step_twists = []
        for step_i in range(n_valid_steps):
            bp_a = step_i                    # first bp of this step
            bp_b = (step_i + 1) % n_bp      # second bp of this step (circular)

            if not (bp_intact[bp_a] and bp_intact[bp_b]):
                continue

            tw = twist[tw_idx, step_i]

            if tw <= 0.0 or tw > 100.0 or tw == 999.0:
                continue

            step_twists.append(tw)

        if step_twists:
            frame_avgs.append(np.mean(step_twists))
            frame_counts.append(len(step_twists))
            included_frames += 1

    frame_avgs   = np.array(frame_avgs)
    frame_counts = np.array(frame_counts)

    if len(frame_avgs) == 0:
        return None

    grand_mean    = np.mean(frame_avgs)
    std_err       = np.std(frame_avgs, ddof=1) / np.sqrt(len(frame_avgs))
    helical_repeat = 360.0 / grand_mean

    return {
        "n_frames_analyzed":   included_frames,
        "mean_twist_deg":      grand_mean,
        "std_err_twist_deg":   std_err,
        "helical_repeat_bpturn": helical_repeat,
        "mean_intact_steps_per_frame": np.mean(frame_counts),
        "frame_avgs": frame_avgs,
    }


# arg parsing

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Compute the average helical repeat of intact base pairs "
            "over the last N (default 10) ns of a circular dsDNA minicircle simulation."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    req = parser.add_argument_group("required arguments")
    req.add_argument("--data_dir",       required=True,
                     help="Directory containing H-bond and stacking data files.")
    req.add_argument("--n_residues",     required=True, type=int,
                     help="Total number of residues (2 * n_bp).")
    req.add_argument("--n_valid_steps",  required=True, type=int,
                     help="Number of valid sequential twist steps in Twist__twist.xvg "
                          "(= n_bp for circular dna).")
    req.add_argument("--dist_xvg",       required=True,
                     help="Name of GROMACS distance .xvg file.")
    req.add_argument("--ang_xvg",        required=True,
                     help="Name of GROMACS angle .xvg file.")
    req.add_argument("--com_dir",        required=True,
                     help="Subdirectory name holding nucleobase COM .xvg files.")
    req.add_argument("--vec_dir",        required=True,
                     help="Subdirectory name holding nucleobase vector .xvg files.")
    req.add_argument("--twist_xvg",      required=True,
                     help="Path to Twist__twist.xvg from do_x3dna.")

    parser.add_argument("--last_ns",          default=10,   type=float,
                        help="Nanoseconds from end of simulation to analyze.")
    parser.add_argument("--hbond_max_dist",   default=0.35, type=float,
                        help="H-bond distance cutoff (nm).")
    parser.add_argument("--hbond_max_angle",  default=30.0, type=float,
                        help="H-bond angle cutoff (degrees).")
    parser.add_argument("--stacking_thresh",  default=0.6,  type=float,
                        help="xi threshold for broken stacking.")

    return parser.parse_args()


# main

def main():
    args = parse_args()
    data_dir = os.path.abspath(args.data_dir)
    print()
    print("  Intact-BP helical repeat analysis")
    print()
    print(f"  Data directory  : {data_dir}")
    print(f"  n_residues      : {args.n_residues}  ({args.n_residues // 2} bp)")
    print(f"  n_valid_steps   : {args.n_valid_steps}")
    print(f"  Last N ns       : {args.last_ns}")
    print(f"  H-bond cutoffs  : dist ≤ {args.hbond_max_dist} nm, "
          f"angle ≤ {args.hbond_max_angle}°")
    print(f"  Stacking cutoff : xi ≤ {args.stacking_thresh}")
    print()

    print("Loading H-bond data …")
    time_hb, hb_intact = get_hbond_intact(
        data_dir, args.dist_xvg, args.ang_xvg,
        args.hbond_max_dist, args.hbond_max_angle,
    )
    print(f"  {len(time_hb)} frames, "
          f"t = {time_hb[0]:.0f} - {time_hb[-1]:.0f} ps")

    print("Loading stacking data ...")
    time_st, st_intact = get_stacking_intact(
        data_dir, args.n_residues, args.com_dir, args.vec_dir,
        args.stacking_thresh,
    )
    print(f"  {len(time_st)} frames")

    print("Loading twist data ...")
    time_tw, twist = load_twist(args.twist_xvg, args.n_valid_steps)
    print(f"  {len(time_tw)} frames, {twist.shape[1]} valid steps loaded")

    print()
    print(f"Analyzing last {args.last_ns} ns ...")
    results = compute_intact_helical_repeat(
        time_hb, hb_intact,
        time_st, st_intact,
        time_tw, twist,
        args.last_ns, args.n_valid_steps,
    )

    if results is None:
        print("ERROR: No frames found in the requested time window, "
              "or no intact steps at any frame.")
        return

    # print diagnostics
    print()
    print()
    print("  Results")
    print()
    print(f"  Frames analyzed           : {results['n_frames_analyzed']}")
    print(f"  Mean intact steps / frame : "
          f"{results['mean_intact_steps_per_frame']:.1f}")
    print(f"  Mean twist (intact steps) : "
          f"{results['mean_twist_deg']:.4f} ± {results['std_err_twist_deg']:.4f} °")
    print(f"  Helical repeat            : "
          f"{results['helical_repeat_bpturn']:.3f} bp/turn")
    print()
    print(f"  Report as: {results['mean_twist_deg']:.3f} ± "
          f"{results['std_err_twist_deg']:.3f} °  "
          f"({results['helical_repeat_bpturn']:.2f} bp/turn)")
    print()


if __name__ == "__main__":
    main()
