# main.jl
# include("XSBench.jl")
include("Simulation.jl")

using Base: write
using .Simulation
import JACC

# using .JACC
# using .XSBench
# using Serialization
# using Profile
# using ProfileView
# Function to check if SD is not null without printing


function main()
    print("XSBench Julia Implementation\n")

    println("Loading simulation data...")
    in = Simulation.LoadData.default_input()
    SD = Simulation.LoadData.load_sim_data()


    println("Calculating macro cross section...")
    Simulation.run_event_based_simulation(in, SD)
    # Simulation.simulation_ex()


end

main()
