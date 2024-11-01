#!/bin/bash

# Set up project directory
PROJ_DIR=/ccs/proj/csc594/ahuante/JACC_APPS/miniBUDEJACC

# good practice to avoid conflicts with existing default modules
module purge

# Load required modules
module load Core/24.07 #needed to load julia in frontier.
module load julia/1.10.2

# Install and set up packages
julia --project=$PROJ_DIR -e 'import Pkg; Pkg.instantiate()'

# Set backend preference
julia --project=$PROJ_DIR -e '
include("src/BasicBUDEPreferences.jl")
using .BasicBUDEPreferences
BasicBUDEPreferences.set_backend("amdgpu")
println("Backend set to: ", BasicBUDEPreferences.get_backend())
'

echo "Configuration complete."

# launch this file with batch `$batch config_frontier_JACC.sh`