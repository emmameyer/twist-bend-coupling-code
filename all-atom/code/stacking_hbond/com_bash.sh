#!/bin/bash

# Create output directory
mkdir -p com_files

# Loop through all xx residues
for i in {1..180}
do
    echo "RESI_${i}" | gmx traj -f  -s  -n nucleobase_COM_atoms.ndx -ox com_files/nucleobase_COM_coord_${i}.xvg -com -ng 1
done
