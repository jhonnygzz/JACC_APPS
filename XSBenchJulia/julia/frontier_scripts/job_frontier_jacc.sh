#!/bin/bash
#SBATCH -A CSC594
#SBATCH -J AMD_XSBench_JACC
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err
#SBATCH -t 00:08:00
#SBATCH -p batch
#SBATCH -N 1

date

XSBench_DIR=/ccs/proj/csc594/jvalglz/JACC_APPS/XSBenchJulia/julia
XSBenchJulia_EXE=$XSBench_DIR/src/main.jl

srun -n 1 --gpus=1 julia --project=$XSBench_DIR $XSBenchJulia_EXE

echo "Job completed."

# launch this file with sbatch `$ sbatch job_frontier_JACC.sh`
