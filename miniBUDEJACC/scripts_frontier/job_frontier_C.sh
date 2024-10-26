#!/bin/bash
#SBATCH -A CSC594
#SBATCH -J AMD_miniBUDE_C
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH -t 00:04:00
#SBATCH -p batch
#SBATCH -N 1

date

# project directory with the executable
cd /ccs/proj/csc594/ahuante/miniBUDE/build

srun -n 1 --gpus=1 ./hip-bude --ppwi 4 --wgsize 64

echo "Job completed."

# launch this file with sbatch `$ sbatch job_frontier_C.sh` AFTER using cmake and building.