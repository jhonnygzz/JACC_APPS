#!/bin/bash

# Set up project directory
PROJ_DIR=/ccs/proj/csc594/ahuante/miniBUDE/src/julia/miniBUDE.jl
KA_DIR=$PROJ_DIR/KernelAbstractions


# good practice to avoid conflicts with existing default modules
module purge

# Load required modules
module load Core/24.00 #needed to load julia in frontier.
module load julia/1.8.2
module load rocm/6.2.0

# Install and set up packages
julia --project=$KA_DIR -e 'import Pkg; Pkg.instantiate()'

echo "Configuration complete."

# launch this file with source `$source config_frontier_JACC.sh`