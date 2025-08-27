#!/bin/bash
#SBATCH -n 1
#SBATCH --gres=gpu:1 						
#SBATCH -o slurm-%A.%a.out 

source /home/software/amber22/amber.sh

echo "CUDA DEVICES:" $CUDA_VISIBLE_DEVICES

srun pmemd.cuda -O -i min.in  -o min.out  -p   system.wat.leap.prmtop  -c system.wat.leap.rst7 -r min.rst7
srun pmemd.cuda -O -i equilibrate1.in  -o eq1.out -p system.wat.leap.prmtop  -c min.rst7  -r eq1.rst7 -x eq1.nc
srun sander -O -i equilibrate2.in   -o eq2.out  -p system.wat.leap.prmtop  -c eq1.rst7  -r eq2.rst7 -x eq2.nc
srun pmemd.cuda -O -i equilibrate2_ext.in -o eq2_ext.out -p system.wat.leap.prmtop -c eq2.rst7 -r eq2_ext.rst7 -x eq2_ext.nc
