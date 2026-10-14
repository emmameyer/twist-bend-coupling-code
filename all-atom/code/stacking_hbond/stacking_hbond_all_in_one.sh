#!/bin/bash
## Overall analysis maker
# define input files
input_gro=*.gro
input_traj=*.xtc
input_tpr=*.tpr
output_location=/Users/your/output/directory
nucleotides=170
bps=85
title='Title'

# Create output directory if it doesn't exist
mkdir -p "$output_location/com_files"
mkdir -p "$output_location/vec_files"

# Change to output directory
cd "$output_location" || exit

# Hbond
printf "Hbond analysis...\n"
printf "Pairs ndx:\n"
python3 /Users/emma/Documents/kent/fraying_analysis/circular_analysis/make_wc_pairs_ndx_circular.py \
"$input_gro" \
"$bps" \
"$output_location"

printf "Triplets ndx:\n"
python3 /Users/emma/Documents/kent/fraying_analysis/circular_analysis/make_wc_triplets_ndx_circular.py \
"$input_gro" \
"$bps" \
"$output_location"

printf "Gmx distance:\n"
gmx distance -f "$input_traj" \
-s "$input_tpr" \
-n wc_pairs.ndx \
-oall hbond_dist.xvg <<EOF
$bps
EOF

printf "Gmx angle:\n"
gmx angle -f "$input_traj" \
-n wc_triplets.ndx \
-type angle -all \
-ov hbond_angle.xvg <<EOF
$bps
EOF

printf "Plotting Hbond analysis...\n"
python3 /Users/emma/Documents/kent/fraying_analysis/circular_analysis/plot_hbond.py \
"HBOND: $title" \
hbond_dist.xvg \
hbond_angle.xvg \
0 \
9.1 \
1.7 \
0 \
"$output_location"

# Stacking
printf "Stacking analysis...\n"
printf "Ndx maker...\n"
python3 /Users/emma/Documents/kent/fraying_analysis/circular_analysis/make_nucleobase_plane_COM_index_files.py \
"$nucleotides" \
"$input_gro" \
"$output_location"

printf "Gmx traj:\n"
printf "Calculating COMs...\n"
for i in $(seq 1 "$nucleotides")
do
echo "RESI_${i}" | gmx traj -f "$input_traj" \
-s "$input_tpr" \
-n nucleobase_COM_atoms.ndx \
-ox com_files/nucleobase_COM_coord_${i}.xvg \
-com -ng 1
done

printf "Calculating vectors... \n"
for i in $(seq 1 "$nucleotides")
do
echo "RESI_${i}" | gmx traj -f "$input_traj" \
-s "$input_tpr" \
-n nucleobase_vec_atoms.ndx \
-ox vec_files/nucleobase_vec_coord_${i}.xvg \
-ng 1
done

printf "Plotting stacking analysis...\n"
python3 /Users/emma/Documents/kent/fraying_analysis/circular_analysis/plot_stacking_colormap.py \
"Stacking: $title" \
com_files \
vec_files \
"$nucleotides" \
1 \
9.1 \
1.7 \
0 \
circular \
"$output_location"

printf "\nAnalysis complete! All files saved to: $output_location\n"
