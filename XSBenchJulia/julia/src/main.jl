
include("Simulation.jl")

using Base: write
using .Simulation


function main()
    print("XSBench Julia Implementation\n")

    println("Loading simulation data...")
    in = Simulation.LoadData.default_input()
    SD = Simulation.LoadData.load_sim_data()
    # println("Calculating macro cross section...")
    Simulation.run_event_based_simulation(in, SD)


    # Simulation.test_simulation()



end


function compare_files(file1_path, file2_path)
    open(file1_path, "r") do file1
        open(file2_path, "r") do file2
            for (line_number, (line1, line2)) in enumerate(zip(eachline(file1), eachline(file2)))
                if line1 != line2
                    println("Difference found at line $line_number:")
                    println("  File cuda/verification_baseline: $line1")
                    println("  File julia/verification_h_data_tests: $line2")
                    return
                end
            end
            
            # Check if one file has more lines than the other
            if !eof(file1)
                println("File cuda/verification_baseline has more lines than File julia/verification_h_data_tests")
            elseif !eof(file2)
                println("File julia/verification_h_data_tests has more lines than File File cuda/verification_baseline")
            else
                println("Files are identical")
            end
        end
    end
end

function count_differences(file1_path, file2_path)
    differences = 0
    open(file1_path, "r") do file1
        open(file2_path, "r") do file2
            for (line1, line2) in zip(eachline(file1), eachline(file2))
                if line1 != line2
                    differences += 1
                end
            end
            
            # Check if one file has more lines than the other
            while !eof(file1)
                readline(file1)
                differences += 1
            end
            while !eof(file2)
                readline(file2)
                differences += 1
            end
        end
    end
    println("Total differences: $differences")
end

main()
compare_files("/home/jvalglz/Fall2024/OLD/XSBench/cuda/verification_baseline_tests.txt", "/home/jvalglz/Fall2024/JACC_APPS/XSBenchJulia/julia/verification_h_data_tests.txt")
count_differences("/home/jvalglz/Fall2024/OLD/XSBench/cuda/verification_baseline_tests.txt", "/home/jvalglz/Fall2024/JACC_APPS/XSBenchJulia/julia/verification_h_data_tests.txt")

