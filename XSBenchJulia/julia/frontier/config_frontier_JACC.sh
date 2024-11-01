#!/bin/bash


# Set up project directory
PROJ_DIR=/ccs/proj/csc594/jvalglz/JACC_APPS/XSBenchJulia
export JULIA_DEPOT_PATH=$PROJ_DIR/.julia
XS_DIR=$PROJ_DIR/julia


# Set up project directory
# PROJ_DIR=/ccs/proj/csc594/jvalglz/JACC_APPS/XSBenchJulia/julia

# remove existing generated Manifest.toml
rm -f $XS_DIR/Manifest.toml
rm -f $XS_DIR/LocalPreferences.toml

# Good practice to avoid conflicts with existing default modules
module purge

# Load required modules
# module load Core/24.07 # Needed to load Julia in Frontier.
# module load rocm/6.2.0
module load PrgEnv-gnu-amd/8.4.0
#module load julia/1.10.4

# doWNLOAD julia from X86, its an executible
# export PATH=$PATH:/ccs/proj/csc594/jvalglz/julia_installs/julia-1.11.1/bin
export PATH=/ccs/proj/csc594/jvalglz/julia_installs/julia-1.11.1/bin:$PATH

# export ROCM_PATH=/opt/rocm-6.2.0/


# Instantiate the project by installing packages in Project.toml
julia --project=$XS_DIR -e 'using Pkg; Pkg.instantiate()'
# Verify the packages are installed correctly
julia --project=$XS_DIR -e 'using Pkg; Pkg.build()'
# JACC.jl and main.jl won't precompile, but other packages will
julia --project=$XS_DIR -e 'using Pkg; Pkg.precompile()'

# Install and set up packages
# julia --project=$PROJ_DIR -e 'import Pkg; Pkg.instantiate()'
# julia --project=$PROJ_DIR -e 'import Pkg; Pkg.add("Preferences")'
# julia --project=$PROJ_DIR -e 'import Pkg; Pkg.add("Example")'
# julia --project=$PROJ_DIR -e 'import Pkg; Pkg.add("AMDGPU")'
# julia --project=$PROJ_DIR -e 'import Pkg; Pkg.add("CUDA")'
# julia --project=$PROJ_DIR -e 'import JACC.JACCPreferences;'
# julia --project=$PROJ_DIR -e 'import Pkg; Pkg.UPDATED_REGISTRY_THIS_SESSION[] = true'julia --project=$PROJ_DIR -e 'import Pkg; Pkg.update()'

# Set backend preference
#julia --project=$PROJ_DIR -e '
#include("src/BasicXSBenchPreferences.jl")
#using .BasicXSBenchPreferences
#BasicXSBenchPreferences.set_backend("jacc-amdgpu")
#println("Backend set to: ", BasicXSBenchPreferences.get_backend())
#'

echo "Configuration complete."

# Launch this file with source `$ source config_frontier_JACC.sh`