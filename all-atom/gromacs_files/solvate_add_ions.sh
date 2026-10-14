#!/bin/bash

# solvate and add ions for general md run on osc
# ensure correct loaded modules
module reset
module load gcc/12.3.0 mvapich/3.0 cuda/12.4.1
module load gromacs/2024.4
module list

# define input gro as the em_1.gro file from previous step
input_gro=$(ls em_1.gro)
input_top=$(ls system.top)
#output_name=${input_top%.top}
output_name=nick_histone
current_dir=$(pwd)
input_ions_mdp=${current_dir}/inputs/ions.mdp
input_em_mdp=${current_dir}/inputs/em.mdp

gmx solvate -cp "$input_gro" -cs spc216.gro -o "${output_name}_solvated.gro" -p "$input_top" 
gmx grompp -f "$input_ions_mdp" -c "${output_name}_solvated.gro" -p "$input_top" -o ions.tpr
gmx genion -s ions.tpr -o "${output_name}_solvated_ions.gro" -p "$input_top" -pname NA -nname CL -neutral -conc 0.150 <<EOF
14
EOF

# second energy minimization after solvation and ion addition
gmx_d grompp -f "$input_em_mdp" -c "${output_name}_solvated_ions.gro" -p "$input_top" -o em_2.tpr
gmx_d mdrun -v -deffnm em_2

# completed check
printf "Solvation and ion addition completed. Output files:\n"
printf "%s\n" "${output_name}_solvated.gro" "${output_name}_solvated_ions.gro" "em_2.tpr" "em_2.gro"
