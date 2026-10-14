#!/bin/bash

# nvt and npt equilibration for general md run on osc
# ensure correct loaded modules
module reset
module load gcc/12.3.0 mvapich/3.0 cuda/12.4.1
module load gromacs/2024.4
module list

# define input gro as the em_2.gro file from previous step
input_gro=$(ls em_2.gro)
input_top=$(ls system.top)
#output_name=${input_top%.top}
output_name=nick_histone
current_dir=$(pwd)
input_nvt_mdp=${current_dir}/inputs/nvt.mdp
input_npt_mdp=${current_dir}/inputs/npt.mdp

gmx grompp -f "$input_nvt_mdp" -c "$input_gro" -r "$input_gro" -p "$input_top" -o nvt.tpr
gmx mdrun -v -deffnm nvt

input_nvt_gro=$(ls nvt.gro)

gmx grompp -f "$input_npt_mdp" -c "$input_nvt_gro" -r "$input_nvt_gro" -p "$input_top" -o npt.tpr
gmx mdrun -v -deffnm npt

printf "NVT and NPT equilibration completed. Output files:\n"
printf "%s\n" "nvt.tpr" "nvt.gro" "npt.tpr" "npt.gro"
