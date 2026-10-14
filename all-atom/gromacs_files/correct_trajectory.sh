#!/bin/bash

# correct trajectories

# input files
input_traj="md_1150ns.xtc"
input_tpr="md_50ns_to_1150ns.tpr"
input_gro="md_50ns_to_1150ns.gro"
input_index="DNA_histone.ndx"
# output name is the same as the .pdb name used in initial setup but without .pdb
#pdb_file=$(ls *.pdb)
#output_name=${pdb_file%.pdb}
output_name=nicked_80bp_minicircle_histone_md_1150ns

gmx trjconv -f $input_traj -s $input_tpr -pbc whole -n $input_ndx.ndx -o step1.xtc << EOF
20
EOF

gmx trjconv -f step1.xtc -s $input_tpr -pbc cluster -n $input_ndx.ndx -o step2.xtc << EOF
20
20
EOF

gmx trjconv -f step2.xtc -s $input_tpr -pbc mol -center -n $input_ndx.ndx -o step3.xtc << EOF
20
20
EOF

gmx trjconv -f step3.xtc -s $input_tpr -fit rot+trans -n $input_ndx.ndx -o ${output_name}.xtc << EOF
20
20
EOF

gmx trjconv -f $input_gro -s $input_tpr -pbc whole -n $input_ndx.ndx -o step1.gro << EOF
20
EOF

gmx trjconv -f step1.gro -s $input_tpr -pbc cluster -n $input_ndx.ndx -o step2.gro << EOF
20
20
EOF

gmx trjconv -f step2.gro -s $input_tpr -pbc mol -center -n $input_ndx.ndx -o step3.gro << EOF
20
20
EOF

gmx trjconv -f step3.gro -s $input_tpr -fit rot+trans -n $input_ndx.ndx -o ${output_name}.gro << EOF
20
20
EOF


# change name of tpr file to match output xtc and gro files, simply for organization purposes
cp $input_tpr ${output_name}.tpr
