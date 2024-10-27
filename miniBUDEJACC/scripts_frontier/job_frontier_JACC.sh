#!/bin/bash
#SBATCH -A CSC594
#SBATCH -J AMD_miniBUDE_JACC
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH -t 00:05:00
#SBATCH -p batch
#SBATCH -N 1

date

BUDE_DIR=/ccs/proj/csc594/ahuante/JACC_APPS/miniBUDEJACC
BUDE_EXE=$BUDE_DIR/src/JACCBUDE.jl

srun -n 1 --gpus=1 julia --project=$BUDE_DIR $BUDE_EXE -p 4

echo "Job completed."

# launch this file with sbatch `$ sbatch job_frontier_JACC.sh`