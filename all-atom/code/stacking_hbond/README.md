# README file for stacking_hbond

## Scripts for building Watson-Crick h-bond index groups and running the full h-bond/stacking interaction analysis for a DNA trajectory.

1. [make_wc_pairs_ndx.py](make_wc_pairs_ndx.py): builds a GROMACS index file with Watson-Crick h-bond pairs (H, A) for gmx distance and triplets (H, D, A) for gmx angle, from a structure file and a pairs.tsv listing which DNA positions pair with which (works for any pairing scheme, not just circular).
2. [make_wc_triplets_ndx.py](make_wc_triplets_ndx.py): same idea as above but only outputs the H-D-A triplet groups for gmx angle.
3. [make_nucleobase_plane_COM_index_files.py](make_nucleobase_plane_COM_index_files.py): creates index files 'nucleobase_vec_atoms.ndx' and 'nucleobase_COM_atoms.ndx' which are used as input for GROMACS utility 'traj'.
4. [com_bash.sh](com_bash.sh): loops over residues and runs gmx traj to get the nucleobase center-of-mass coordinates for each, used in stacking analysis. Needs a nucleobase_COM_atoms.ndx file and the residue range edited per system, using Rachel Bricker's script to create it.
5. [vec_bash.sh](vec_bash.sh): same as com_bash.sh but pulls the nucleobase plane vectors instead of COM, needs nucleobase_vec_atoms.ndx.
6. [stacking_hbond_all_in_one.sh](stacking_hbond_all_in_one.sh): chains the full pipeline together for a circular dsDNA system, builds the pairs/triplets ndx files, runs gmx distance and gmx angle for h-bonds, plots the h-bond results, builds the nucleobase COM index, calculates COM and vector coordinates per residue, and plots the stacking colormap. Edit the input file paths, output_location, nucleotides, bps, and title at the top before running, and note the script paths inside point to a local fraying_analysis directory that will need updating to match wherever those scripts actually live.

## Hydrogen bond analysis
1. create the index files for the pairs and triplets using the python scripts.
2. run gmx angle using the triplet index, outputs one xvg file (bash shortcut found in [the overall analysis script](stacking_hbond_all_in_one.sh)).
    * gmx angle -f input.xtc -n wc_triplets.ndx -type angle -all -ov hbond_angle.xvg
3. run gmx distance using the pair index, outputs one xvg file (bash shortcut found in [the overall analysis script](stacking_hbond_all_in_one.sh)).
    * gmx distance -f input.xtc -s input.tpr -n wc_pairs.ndx -oall hbond_dist.xvg
4. analyze and plot using [plot_hbond.py](../plotters/plot_hbond.py)

## Stacking interaction analysis
1. create index files using Rachel Bricker's make_nucleobase_plane_COM_index_files.py
2. run gmx traj for COM and vec ndx files, use the bash script shortcuts to run through all nucleotides, found in [com_bash](com_bash), [vec_bash](vec_bash), or [the overall analysis script](stacking_hbond_all_in_one.sh).
3. analyze and plot using [plot_stacking.py](../plotters/plot_stacking.py)
