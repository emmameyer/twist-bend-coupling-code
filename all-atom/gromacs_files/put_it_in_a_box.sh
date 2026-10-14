#!/bin/bash
# initial steps for general md run on osc

# ensure correct loaded modules
module reset
module load gcc/12.3.0 mvapich/3.0 cuda/12.4.1
module load gromacs/2024.4
module list

# define pdb as only pdb in directory
#input_pdb=$(ls *.pdb)
input_ff=$(amber14sb_cufix_OL21)

# define output file names to be the same as input pdb but without .pdb extension
#output_name=${input_pdb%.pdb}
output_name1=system
output_name=nicked_histone

# ensure output directorie is same as input directory
output_location=$(pwd)

gmx editconf -f "${output_name1}_processed.gro" -o "${output_name}_newbox.gro" -c -d 1.0 -bt cubic

# completed
printf "It's now in a box!"
