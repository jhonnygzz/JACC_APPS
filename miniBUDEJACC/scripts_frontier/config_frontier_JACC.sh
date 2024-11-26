#!/bin/bash

# Set up project directory
PROJ_DIR=/ccs/proj/csc594/ahuante/JACC_APPS/miniBUDEJACC

# good practice to avoid conflicts with existing default modules
module purge

# Load required modules
module load Core/24.07 #needed to load julia in frontier.
module load julia/1.10.2

# Set backend preference
julia --project=$PROJ_DIR -e '
using JACC.JACCPreferences
JACCPreferences.set_backend("amdgpu")
using JACC
@show JACC.JACCPreferences.backend  # Should show "amdgpu"
'

# Install and set up packages
julia --project=$PROJ_DIR -e 'import Pkg; Pkg.instantiate()'

echo "Configuration complete."

# launch this file with source `$source config_frontier_JACC.sh`