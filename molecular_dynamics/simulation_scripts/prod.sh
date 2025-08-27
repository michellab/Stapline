#!/bin/bash
#SBATCH -n 1
#SBATCH --gres=gpu:1 						
#SBATCH -o slurm-%A.%a.out 

source /home/software/amber22/amber.sh

echo "CUDA DEVICES:" $CUDA_VISIBLE_DEVICES

srun pmemd.cuda -O -i prod.in -o prod.out -x prod.nc -p system.wat.leap.prmtop -c eq2_ext.rst7 -r prod.rst7