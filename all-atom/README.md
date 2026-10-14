# All-atom MD simulation and analysis code

GROMACS setup inputs and Python/shell analysis scripts for the all-atom
simulations in "Twist-Dependent Elastic Limits and Unwinding of Tightly
Bent DNA" (*Nucleic Acids Research*, 2026). Authors, citation, and DOI:
see the top-level README.

## Software

| Software | Version |
|---|---|
| GROMACS (simulations) | 2024.4 |
| do_x3dna | 2025.2-dev (installed alongside GROMACS 2025) |
| 3DNA (used by do_x3dna) |  2.4.8-2023 |
| Python | 3.12.11 |
| NumPy | 2.3.1 |
| SciPy | 1.15.2 |
| Matplotlib | 3.10.3 |
| MDAnalysis | 2.9.0 |
| PyMOL | 3.1.0 |

Force field: `amber14sb_parmbsc1_cufix` with TIP3P water, 0.150 M NaCl.
Production simulations were run on the Ohio Supercomputer Center.

## Contents

- [gromacs_files](gromacs_files): shell scripts and `.mdp` files for the MD pipeline
  (`pdb2gmx` through production run). See [gromacs_files/README.md](gromacs_files/README.md).
- [do_x3dna](do_x3dna): the do_x3dna command used to extract base-pair step
  parameters, including twist.
- [code/twist_analysis](code/twist_analysis): extracts and block-averages twist angles from
  do_x3dna output.
- [code/stacking_hbond](code/stacking_hbond): builds Watson-Crick h-bond index groups and
  analyzes h-bonding and base stacking.
- [code/plotters](code/plotters): scripts for the figures (twist, h-bonds, stacking,
  free energy of circularization, helical repeat) plus shared styling
  (`fig_style.py`).
- [code/intact_bps_h_analysis.py](code/intact_bps_h_analysis.py): average helical repeat of intact base
  pairs over the last N ns.
- [code/merge_gro_files.py](code/merge_gro_files.py): merges `.gro` files (e.g. DNA + histone).

## Workflow

1. Build and run the system with [gromacs_files](gromacs_files) (setup parameters below).
2. Extract base-pair step parameters with the command in
   [do_x3dna/twist_analysis.txt](do_x3dna/twist_analysis.txt). The key output is the `Twist_*.xvg` file.
3. Analyze:
   - Twist: run [code/twist_analysis/process_twist.py](code/twist_analysis/process_twist.py) on `Twist_*.xvg`,
     then plot with [code/plotters/plot_h.py](code/plotters/plot_h.py) (or block-average with
     `block_average.py` and plot with `plot_twist.py`).
   - H-bonds and stacking: follow [code/stacking_hbond/README.md](code/stacking_hbond/README.md).
4. Plot with [code/plotters](code/plotters) (see [code/plotters/README.md](code/plotters/README.md)).

Many scripts are hard-coded for a given system (file names, residue
ranges, number of base pairs) and need editing before use. These are
noted in each script and README.

## Simulation setup (summary)

Energy minimization (double precision), solvation in a cubic box
(1.0 nm padding), neutralization and 0.150 M NaCl, second minimization,
NVT then NPT equilibration, then production. Parameters are in
`gromacs_files/inputs/*.mdp`.

## Notes

- [code/plotters/plot_sawtooths.py](code/plotters/plot_sawtooths.py) and related scripts read `.txt` files
  digitized from the gel data shown in the paper. These input files are
  not included.
- The stacking analysis uses `make_nucleobase_plane_COM_index_files.py`
  (R. Bricker), included in [code/stacking_hbond](code/stacking_hbond/).
- Trajectories are not included in this archive.

## Credits

- do_x3dna: R. Kumar, https://github.com/rjdkmr/do_x3dna
- 3DNA: Lu & Olson, *Nucleic Acids Res.* 31, 5108 (2003)
- GROMACS: https://www.gromacs.org

## License

GNU GPL v3.0 (see the top-level `LICENSE`).
