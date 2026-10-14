#!/bin/bash

# Create output directory
mkdir -p vec_files

# Loop through all xx residues
for i in {1..180}
do
    echo "RESI_${i}" | gmx traj -f *.xtc -s *.tpr -n nucleobase_vec_atoms.ndx -ox vec_files/nucleobase_vec_coord_${i}.xvg -ng 1
done
