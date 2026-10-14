#!/bin/bash

# first energy minimization for general md run on osc
# ensure correct loaded modules
module reset
module load gcc/12.3.0 mvapich/3.0 cuda/12.4.1
module load gromacs/2024.4
module list

# define input gro as only gro in directory ending with _newbox.gro
input_gro=$(ls *_newbox.gro)
input_top=$(ls system.top)
#output_name=${input_top%.top}
output_name=nick_histone
current_dir=$(pwd)
input_em_mdp=${current_dir}/inputs/em.mdp

gmx_d grompp -f "$input_em_mdp" -c "$input_gro" -p "$input_top" -o em_1.tpr -maxwarn 1
# run energy minimization
gmx_d mdrun -v -deffnm em_1

# completed
printf "First energy minimization completed!"
