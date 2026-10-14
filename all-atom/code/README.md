README file for code

This folder holds the analysis pipeline: extracting twist, h-bond, and stacking data from trajectories, finding which base pairs stay intact, and plotting results.

Subfolders:
1. [twist_analysis](twist_analysis): extracting and block-averaging twist angles from do_x3dna output.
2. [stacking_hbond](stacking_hbond): builds Watson-Crick h-bond index groups (pairs/triplets) and runs the COM/vector/h-bond/stacking analysis for a trajectory.
3. [plotters](plotters): plotting scripts for twist, h-bonds, stacking, free energy of circularization, and shared figure styling.

Standalone scripts:
1. [intact_bps_h_analysis.py](intact_bps_h_analysis.py): finds the average helical repeat (bp/turn) of only the intact base pairs of a circular dsDNA over the last N ns of a simulation, where a base pair is intact if both the h-bond (dist <= 0.35 nm, angle <= 30 deg) and stacking (xi <= 0.6) criteria are met for both residues. See the script's docstring for full usage and arguments.
2. [merge_gro_files.py](merge_gro_files.py): merges two .gro files (e.g. DNA + histone) into one combined system using MDAnalysis.
