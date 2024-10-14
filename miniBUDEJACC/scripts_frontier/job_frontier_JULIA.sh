#!/bin/bash
#SBATCH -A CSC594
#SBATCH -J AMD_miniBUDE_JULIA
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH -t 00:04:00
#SBATCH -p batch
#SBATCH -N 1

date

BUDE_DIR=/ccs/proj/csc594/ahuante/miniBUDE/src/julia/miniBUDE.jl
BUDE_AMDGPU=$BUDE_DIR/AMDGPU
BUDE_EXE=$BUDE_DIR/src/AMDGPU.jl

srun -n 1 --gpus=1 julia --project=$BUDE_AMDGPU $BUDE_EXE

echo "Job completed."

# launch this file with sbatch `$ sbatch job_frontier_JACC.sh`