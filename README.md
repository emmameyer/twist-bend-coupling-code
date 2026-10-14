# Code for: Twist-Dependent Elastic Limits and Unwinding of Tightly Bent DNA

Analysis code accompanying the paper "Twist-Dependent Elastic Limits and
Unwinding of Tightly Bent DNA" (*Nucleic Acids Research*, 2026).

Archived at Zenodo: [https://doi.org/10.5281/zenodo.XXXXXXX CHANGE]

## Authors

Soumya Chandrasekhar, Emma Meyer, Fatemeh Fadaei,
Rachel Bricker, Thomas P. Swope, Morgan Schreck, Mark Macri,
Draven Houser, Daniel R. Hollis, John J. Portman, Thorsten L. Schmidt

## Contents

- `oxdna/`: analysis scripts for the coarse-grained oxDNA simulations
  (kink detection, bend and twist angles, kymographs). See
  `oxdna/README.md`.
- `all-atom/`: GROMACS setup inputs and analysis scripts for the
  all-atom simulations (twist, hydrogen bonding, base stacking). See
  `all-atom/README.md`.

This repository contains analysis code and simulation setup files only.
Simulation trajectories are not included.

## Software

Software requirements and credits for third-party tools are listed in
the README of each folder. The oxDNA analysis builds on the oxDNA
analysis tools (OAT) (https://github.com/lorenzo-rovigatti/oxDNA); the
all-atom analysis uses GROMACS, do_x3dna, and MDAnalysis.

## How to cite

If you use this code, please cite the paper and this archive
(see `CITATION.cff` or the "Cite this repository" button on GitHub).

## License

GNU General Public License v3.0 (see `LICENSE`). Some oxDNA files are
modified from the GPL-licensed oxDNA analysis tools; modifications are
noted.
