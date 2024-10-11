#!/bin/bash

# Set up project directory
PROJ_DIR=$HOME/Fall2024/JACC_APPS/miniBUDEJACC
#GS_DIR=$PROJ_DIR/miniBUDEJACC

# remove existing generated Manifest.toml
rm -f $PROJ_DIR/Manifest.toml
rm -f $PROJ_DIR/LocalPreferences.toml

# good practice to avoid conflicts with existing default modules
module purge

# Load required modules
module load julia/1.10.4
module load rocm/6.0.2

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
