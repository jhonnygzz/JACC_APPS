#!/bin/bash
#SBATCH -A CSC594
#SBATCH -J AMD_miniBUDE_C
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH -t 00:05:00
#SBATCH -p batch
#SBATCH -N 1

date

# project directory with the executable
cd /ccs/proj/csc594/ahuante/miniBUDE/build_amdMI250X

srun -n 1 --gpus=1 ./hip-bude --deck ../data/bm2 --ppwi 1 --wgsize 128

echo "Job completed."

# launch this file with sbatch `$ sbatch job_frontier_C.sh` AFTER using cmake and building.