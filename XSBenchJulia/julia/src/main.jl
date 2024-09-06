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
    # in = Simulation.LoadData.default_input()
    # SD = Simulation.LoadData.load_sim_data()

    # Simulation.hello(in)



    println("Calculating macro cross section...")
    # Simulation.run_event_based_simulation(in, SD)
    Simulation.simulation_ex()


end

function mainold()
    # Profile.clear()
    # @profile begin
    # @time begin
    p_energy = 0.626011
    mat = 5
    n_iso = 355
    in_grid = 11303
    in_gridtype = 0
    in_hashbins = 10000
    macro_xs_vector = zeros(Float64, 355)

    println("XSBench Simulation")
    SD = Simulation.init_simulation_data()

    println("Calculating macro cross section...")
    Simulation.calculate_macro_xs(p_energy, mat, n_iso, in_grid, SD.num_nucs, SD.concs, SD.unionized_energy_array, SD.index_grid, SD.nuclide_grid, SD.mats, macro_xs_vector, in_gridtype, in_hashbins, SD.max_num_nucs)
    check_SD_data(SD)
    
    # Profile.print(format=:flat, sortedby=count)
    # ProfileView.view()

    end



function axpy(i, alpha, x, y)
    if i <= length(x)
    @inbounds x[i] += alpha * y[i]
    end
end

function JACCTEST()
    print("JACC Test\n")
    N = 100000000
    # Generate random vectors x and y of length N for the interval [0, 100]
    x = round.(rand(Float32, N) * 100)
    y = round.(rand(Float32, N) * 100)
    alpha = 2.5

    x_d = JACC.Array(x)
    y_d = JACC.Array(y)
    JACC.parallel_for(N, axpy, alpha, x_d, y_d)
    print("Done\n")
    
end

function JACCXSBench()
    p_energy = 0.626011
    mat = 5
    n_iso = 355
    in_grid = 11303
    in_gridtype = 0
    in_hashbins = 10000
    macro_xs_vector = zeros(Float64, 355)

    println("XSBench Simulation")
    SD = Simulation.init_simulation_data()

    println("Calculating macro cross section...")
    # Simulation.calculate_macro_xs(p_energy, mat, n_iso, in_grid, SD.num_nucs, SD.concs, SD.unionized_energy_array, SD.index_grid, SD.nuclide_grid, SD.mats, macro_xs_vector, in_gridtype, in_hashbins, SD.max_num_nucs)
    JACC.parallel_for(1, Simulation.calculate_macro_xs, p_energy, mat, n_iso, in_grid, SD.num_nucs, SD.concs, SD.unionized_energy_array, SD.index_grid, SD.nuclide_grid, SD.mats, macro_xs_vector, in_gridtype, in_hashbins, SD.max_num_nucs)

    check_SD_data(SD)
end



main()
