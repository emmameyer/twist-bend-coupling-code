#!/bin/bash
#SBATCH --job-name=semicircle_histone   # Job name specific to this run
#SBATCH -o %x.o%j
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=emeyer12@kent.edu
#SBATCH --time=160:00:00
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=12 
#SBATCH --gpus-per-node=1
#SBATCH --account=PGS0408

set +vx
module reset
module load gcc/12.3.0 mvapich/3.0 cuda/12.4.1
module load gromacs/2024.4
module list
set -vx

nvidia-smi

# Move to the job submission directory
cd $SLURM_SUBMIT_DIR

# ensure needed files are inside the temporary directory
cp -r inputs $TMPDIR/
cp -r amber14sb_cufix_OL21.ff $TMPDIR/
cp -r system.top posre.itp em* npt* nvt* $TMPDIR/ || cp -r * $TMPDIR/

# Move into the temporary directory for the run
cd $TMPDIR
pwd
ls -l

# continue from npt step
gmx_gpu grompp -f inputs/md.mdp -c npt.gro -t npt.cpt -r npt.gro -p system.top -o md_100ns.tpr

gmx_gpu mdrun -v -deffnm md_100ns -ntomp 12 -nb gpu -pme gpu -bonded gpu

# comment out the above lines and uncomment the below to continue from a previous production run
#gmx_gpu grompp -f inputs/md.mdp -c md_100ns.gro -t md_100ns.cpt -r md_100ns.gro -p system.top -o md_100_to_200ns.tpr

#gmx_gpu mdrun -v -deffnm md_100_to_200ns -ntomp 12 -nb gpu -pme gpu -bonded gpu

# copy the contents of this temp directory to the actual directory
cp -p * $SLURM_SUBMIT_DIR/
ls -al
