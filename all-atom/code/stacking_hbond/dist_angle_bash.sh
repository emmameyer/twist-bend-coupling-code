#!/bin/bash

input_xtc=$(ls *.xtc)
input_tpr=$(ls *.tpr)


gmx distance -s $input_tpr -f $input_xtc -n wc_pairs.ndx -oall hbond_dist.xvg << EOF
all_pairs
EOF

gmx angle -f $input_xtc -n wc_triplets.ndx -type angle -all -ov hbond_angle.xvg << EOF
all_triplets
EOF
