#!/bin/bash

# Set up project directory
PROJ_DIR=/ccs/proj/csc594/jvalglz/JACC_APPS/XSBenchJulia/julia

# Good practice to avoid conflicts with existing default modules
module purge

# Load required modules
module load Core/24.07 # Needed to load Julia in Frontier.
module load julia/1.10.2
module load rocm/6.2.0

# Install and set up packages
julia --project=$PROJ_DIR -e 'import Pkg; Pkg.instantiate()'
julia --project=$PROJ_DIR -e 'import Pkg; Pkg.add("Preferences")'
julia --project=$PROJ_DIR -e 'import Pkg; Pkg.add("Example")'
# Set backend preference
julia --project=$PROJ_DIR -e '
include("src/BasicXSBenchPreferences.jl")
using .BasicXSBenchPreferences
BasicXSBenchPreferences.set_backend("jacc_amdgpu")
println("Backend set to: ", BasicXSBenchPreferences.get_backend())
'

echo "Configuration complete."

# Launch this file with source `$ source config_frontier_JACC.sh`
