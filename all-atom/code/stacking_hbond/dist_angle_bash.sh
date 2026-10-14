#!/bin/bash

input_xtc=$(ls /Users/emma/Documents/kent/minicircles/OL21/md_files/G_*.xtc)
input_tpr=$(ls /Users/emma/Documents/kent/minicircles/OL21/md_files/G_*.tpr)


gmx distance -s $input_tpr -f $input_xtc -n wc_pairs.ndx -oall hbond_dist.xvg << EOF
all_pairs
EOF

gmx angle -f $input_xtc -n wc_triplets.ndx -type angle -all -ov hbond_angle.xvg << EOF
all_triplets
EOF
