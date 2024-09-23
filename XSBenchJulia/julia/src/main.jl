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

    # Print all data inside SD

    # print

    
    
    # Simulation.simulation_ex()

    # Simulation.test_simulation()



end



# main()


function compare_files(file1_path, file2_path)
    open(file1_path, "r") do file1
        open(file2_path, "r") do file2
            for (line_number, (line1, line2)) in enumerate(zip(eachline(file1), eachline(file2)))
                if line1 != line2
                    println("Difference found at line $line_number:")
                    println("  File 1: $line1")
                    println("  File 2: $line2")
                    return
                end
            end
            
            # Check if one file has more lines than the other
            if !eof(file1)
                println("File 1 has more lines than File 2")
            elseif !eof(file2)
                println("File 2 has more lines than File 1")
            else
                println("Files are identical")
            end
        end
    end
end


main()
# compare_files("/home/jvalglz/Fall2024/OLD/XSBench/cuda/verification_baseline.txt", "/home/jvalglz/Fall2024/JACC_APPS/XSBenchJulia/julia/verification_h_data.txt")
