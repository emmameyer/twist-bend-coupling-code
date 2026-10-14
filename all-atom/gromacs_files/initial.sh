#!/bin/bash
# initial steps for general md run on osc

# ensure correct loaded modules
module reset
module load gcc/12.3.0 mvapich/3.0 cuda/12.4.1
module load gromacs/2024.4
module list

# define pdb as only pdb in directory
input_dna_pdb=$(ls *1.pdb)
input_dna_pdb2=$(ls *2.pdb)
input_histone_pdb=$(ls histone*.pdb)
input_ff=$(amber14sb_cufix_OL21)

# define output file names to be the same as input pdb but without .pdb extension
output_name1=${input_dna_pdb%.pdb}
output_name2=${input_dna_pdb2%.pdb}
output_name2=${input_histone_pdb%.pdb}

# ensure output directorie is same as input directory
output_location=$(pwd)

gmx pdb2gmx -f "$input_dna_pdb" -o "${output_name1}_processed.gro" -p "${output_name1}.top" -i "${output_name1}_merged.itp" -water tip3p -ff amber14sb_cufix_OL21 -merge all
gmx pdb2gmx -f "$input_dna_pdb2" -o "${output_name2}_processed.gro" -p "${output_name2}.top" -i "${output_name2}_merged.itp" -water tip3p -ff amber14sb_cufix_OL21 -merge all
gmx pdb2gmx -f "$input_histone_pdb" -o "${output_name3}_processed.gro" -p "${output_name3}.top" -i "${output_name3}_merged.itp" -water tip3p -ff amber14sb_cufix_OL21 -merge all

#gmx pdb2gmx -f *.pdb -o *_processed.gro -p *.top -water tip3p -ff amber14sb_cufix_OL21 -merge all

# completed
printf "All done with pdb2gmx!"
