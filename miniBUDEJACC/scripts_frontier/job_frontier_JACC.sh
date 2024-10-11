#!/bin/bash
#SBATCH -A CSC594
#SBATCH -J AMD_miniBUDE_JACC
#SBATCH -o %x-%j.out
#SBATCH -t 00:01:00
#SBATCH -p batch
#SBATCH -N 1

date

GS_DIR=$HOME/Fall2024/JACC_APPS/miniBUDEJACC
GS_EXE=$GS_DIR/src/JACCBUDE.jl

srun -n 1 --gpus=1 julia --project=$GS_DIR $GS_EXE

echo "Job completed."

# launch this file with sbatch `$ sbatch job_frontier_JACC.sh`