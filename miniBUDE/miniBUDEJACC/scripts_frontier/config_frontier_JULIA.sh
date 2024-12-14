#!/bin/bash

# Set up project directory
PROJ_DIR=/ccs/proj/csc594/ahuante/miniBUDE/src/julia/miniBUDE.jl
AMDGPU_DIR=$PROJ_DIR/AMDGPU


# good practice to avoid conflicts with existing default modules
module purge

# Load required modules
module load Core/24.07 #needed to load julia in frontier.
module load julia/1.10.2
module load rocm/6.2.0

# Install and set up packages
julia --project=$AMDGPU_DIR -e 'import Pkg; Pkg.instantiate()'

echo "Configuration complete."

# launch this file with source `$source config_frontier_JACC.sh`