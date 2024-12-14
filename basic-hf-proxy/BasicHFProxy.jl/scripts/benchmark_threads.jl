# benchmark_threads.jl

using BasicHFProxy
using Statistics

# Test parameters matching Fortran version
#const natom = 8192
const natom = 1024
const ngauss = 3
const txpnt = [6.3624214, 1.1589230, 0.3136498]
const tcoef = [0.154328967295, 0.535328142282, 0.444634542185]

# Generate geometry data using same random number generator as Fortran
# function generate_geometry(natom)
#     tgeom = zeros(Float64, natom * 3)
#     seed = 12345
#     for i in 1:natom*3
#         seed = mod(1103515245 * seed + 12345, Int64(2)^31)
#         n = mod(seed, 181)
#         tgeom[i] = (n / 10.0) * BasicHFProxy.tobohrs
#     end
#     return tgeom
# end

#Changed to 3×natom matrix instead of vector
function generate_geometry(natom)
    tgeom = zeros(Float64, 3, natom)
    seed = 12345
    for i in 1:natom
        for j in 1:3
            seed = mod(1103515245 * seed + 12345, Int64(2)^31)
            n = mod(seed, 181)
            tgeom[j, i] = (n / 10.0) * BasicHFProxy.tobohrs
        end
    end
    return tgeom
end
# function generate_geometry(natom)
#     tgeom = zeros(Float64, natom * 3)  # 1D array, just like Fortran
#     seed = 12345
#     for i in 1:(natom*3)
#         seed = mod(1103515245 * seed + 12345, Int64(2)^31)
#         n = mod(seed, 181)
#         tgeom[i] = (n / 10.0) * BasicHFProxy.tobohrs
#     end
#     geom = zeros(Float64, 3, natom)  # 2D array for computation
#     for i in 1:natom
#         geom[1,i] = tgeom[3i-2]
#         geom[2,i] = tgeom[3i-1]
#         geom[3,i] = tgeom[3i]
#     end
#     return geom  # Return the 2D array
# end


# Create input data structure
tgeom = generate_geometry(natom)
input = (; ngauss, natom, xpnt=txpnt, coef=tcoef, geom=tgeom)

# Benchmark function
function run_benchmark(func, input)
    println("Performing warmup run...")
    E, t = func(input)
    println("Warmup E = ", E)
    println("Warmup section time: ", round(t; digits=5), " seconds")
    
    println("\nPerforming 10 timed runs:")
    times = Float64[]
    energies = Float64[]
    sum_total = 0.0
    
    for i in 1:10
        E, t = func(input)

        println("Run ", i, ": Time = ", round(t; digits=5), " seconds, 2e- energy = ", E)
        sum_total += t
        flush(stdout)
    end
    
    avg_time = sum_total / 10
    println("\nAverage time of 10 calls: ", round(avg_time; digits=5), " seconds")
end

# Run benchmarks
printstyled("Running Atomix version:\n"; bold=true, color=:cyan)
run_benchmark(bhfp_threads_atomix, input)

printstyled("\nRunning Lock version:\n"; bold=true, color=:cyan)
run_benchmark(bhfp_threads_lock, input)
