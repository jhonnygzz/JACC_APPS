module LoadData
# import Simulation

export NuclideGridPoint, default_input, Input, immutableSimulationData, immutableNuclideGridPoint


struct Input
    nthreads::Int
    n_isotopes::Int
    n_gridpoints::Int
    lookups::Int
    # HM::String
    # HM::Vector{UInt8}
    HM::Char
    grid_type::Int
    hash_bins::Int
    particles::Int
    simulation_method::Int
    binary_mode::Int
    kernel_id::Int
end



mutable struct NuclideGridPoint
    energy::Float64
    total_xs::Float64
    elastic_xs::Float64
    absorbtion_xs::Float64
    fission_xs::Float64
    nu_fission_xs::Float64
end

mutable struct SimulationData
    num_nucs::Vector{Int64}
    concs::Vector{Float64}
    mats::Vector{Int64}
    unionized_energy_array::Vector{Float64}
    index_grid::Vector{Int64}
    nuclide_grid::Vector{NuclideGridPoint}
    length_num_nucs::Int32
    length_concs::Int32
    length_mats::Int32
    length_unionized_energy_array::Int32
    length_index_grid::Int64
    length_nuclide_grid::Int32
    max_num_nucs::Int64
    p_energy_samples::Vector{Float64}
    length_p_energy_samples::Int32
    mat_samples::Vector{Int32}
    length_mat_samples::Int32
end

struct immutableSimulationDataOld
    num_nucs::Vector{Int64}
    concs::Vector{Float64}
    mats::Vector{Int64}
    unionized_energy_array::Vector{Float64}
    index_grid::Vector{Int64}
    nuclide_grid::Vector{NuclideGridPoint}
    length_num_nucs::Int32
    length_concs::Int32
    length_mats::Int32
    length_unionized_energy_array::Int32
    length_index_grid::Int64
    length_nuclide_grid::Int32
    max_num_nucs::Int64
    p_energy_samples::Vector{Float64}
    length_p_energy_samples::Int32
    mat_samples::Vector{Int32}
    length_mat_samples::Int32
end

struct immutableNuclideGridPoint
    energy::Float64
    total_xs::Float64
    elastic_xs::Float64
    absorbtion_xs::Float64
    fission_xs::Float64
    nu_fission_xs::Float64
end

struct immutableSimulationData
    num_nucs::Vector{Int64}
    concs::Vector{Float64}
    mats::Vector{Int64}
    unionized_energy_array::Vector{Float64}
    index_grid::Vector{Int64}
    nuclide_grid::Vector{immutableNuclideGridPoint}
    length_num_nucs::Int32
    length_concs::Int32
    length_mats::Int32
    length_unionized_energy_array::Int32
    length_index_grid::Int64
    length_nuclide_grid::Int32
    max_num_nucs::Int64
    p_energy_samples::Vector{Float64}
    length_p_energy_samples::Int32
    mat_samples::Vector{Int32}
    length_mat_samples::Int32
end





struct GPUImmutableNuclideGridPoint
    energy::Float64
    total_xs::Float64
    elastic_xs::Float64
    absorbtion_xs::Float64
    fission_xs::Float64
    nu_fission_xs::Float64
end


struct GPUImmutableSimulationData
    num_nucs::Ptr{Int64}
    concs::Ptr{Float64}
    mats::Ptr{Int64}
    unionized_energy_array::Ptr{Float64}
    index_grid::Ptr{Int64}
    nuclide_grid::Ptr{GPUImmutableNuclideGridPoint}
    length_num_nucs::Int32
    length_concs::Int32
    length_mats::Int32
    length_unionized_energy_array::Int32
    length_index_grid::Int64
    length_nuclide_grid::Int32
    max_num_nucs::Int64
    p_energy_samples::Ptr{Float64}
    length_p_energy_samples::Int32
    mat_samples::Ptr{Int32}
    length_mat_samples::Int32
end



function convert_to_gpu_data(data::immutableSimulationData)::GPUImmutableSimulationData
    gpu_nuclide_grid = CUDA.fill(GPUImmutableNuclideGridPoint(0.0, 0.0, 0.0, 0.0, 0.0, 0.0), length(data.nuclide_grid))
    for i in 1:length(data.nuclide_grid)
        gpu_nuclide_grid[i] = convert_to_gpu_nuclide_grid_point(data.nuclide_grid[i])
    end
    return GPUImmutableSimulationData(
        CUDA.pointer(data.num_nucs),
        CUDA.pointer(data.concs),
        CUDA.pointer(data.mats),
        CUDA.pointer(data.unionized_energy_array),
        CUDA.pointer(data.index_grid),
        CUDA.pointer(gpu_nuclide_grid),
        data.length_num_nucs,
        data.length_concs,
        data.length_mats,
        data.length_unionized_energy_array,
        data.length_index_grid,
        data.length_nuclide_grid,
        data.max_num_nucs,
        CUDA.pointer(data.p_energy_samples),
        data.length_p_energy_samples,
        CUDA.pointer(data.mat_samples),
        data.length_mat_samples
    )
end

function convert_to_gpu_nuclide_grid_point(point::immutableNuclideGridPoint)::GPUImmutableNuclideGridPoint
    return GPUImmutableNuclideGridPoint(
        point.energy,
        point.total_xs,
        point.elastic_xs,
        point.absorbtion_xs,
        point.fission_xs,
        point.nu_fission_xs
    )
end


function default_input()
    return Input(
        1,          # nthreads
        355,        # n_isotopes
        11303,      # n_gridpoints
        2,   # lookups
        'l',        # HM originally "large"
        0,          # grid_type
        10000,      # hash_bins
        0,          # particles
        2,          # simulation_method
        0,          # binary_mode
        0           # kernel_id
    )
end

function load_sim_data()
    fname = "../../../OLD/XSBench/cuda/SimulationData.txt"
    println("Attempting to open file $fname...")


    # Attempt to open the file in read mode
    file = open(fname, "r")

    # Check if the file was successfully opened
    if file === nothing
        error("Failed to open file: $fname")
    end

    println("Reading Data from file...")
    num_nucs = Vector{Int32}(undef, 12)
    concs = Vector{Float64}(undef, 3852)
    mats = Vector{Int32}(undef, 3852)
    unionized_energy_array = Vector{Float64}(undef, 4012565)
    index_grid = Vector{Int32}(undef, 1424460575)
    nuclide_grid = Vector{immutableNuclideGridPoint}(undef, 4012565)
    p_energy_samples = Vector{Float64}(undef, 0)
    mat_samples = Vector{Int32}(undef, 0)
    length_num_nucs = 12
    length_concs = 3852
    length_mats = 3852
    length_unionized_energy_array = 4012565
    length_index_grid = 1424460575
    length_nuclide_grid = 4012565
    length_p_energy_samples = 0
    length_mat_samples = 0
    max_num_nucs = 321

    # Read the num_nucs value amount of integers from the file
    for i in 1:length_num_nucs
        num_nucs[i] = parse(Int32, readline(file))
    end

    # Print all data inside num_nucs array
    println("num_nucs: ", num_nucs)


    # Read the concs value amount of floats from the file
    for i in 1:length_concs
        concs[i] = parse(Float64, readline(file))
    end

    # Read the mats value amount of integers from the file
    for i in 1:length_mats
        mats[i] = parse(Int32, readline(file))
    end

    # Read the unionized_energy_array value amount of floats from the file
    for i in 1:length_unionized_energy_array
        unionized_energy_array[i] = parse(Float64, readline(file))
    end

    # Read the index_grid value amount of integers from the file
    for i in 1:length_index_grid
        index_grid[i] = parse(Int32, readline(file))
    end

    # Read the nuclide_grid value amount of NuclideGridPoint structs from the file
    for i in 1:length_nuclide_grid
        line = readline(file)
        energy = parse(Float64, line)
        line = readline(file)
        total_xs = parse(Float64, line)
        line = readline(file)
        elastic_xs = parse(Float64, line)
        line = readline(file)
        absorbtion_xs = parse(Float64, line)
        line = readline(file)
        fission_xs = parse(Float64, line)
        line = readline(file)
        nu_fission_xs = parse(Float64, line)
        nuclide_grid[i] = immutableNuclideGridPoint(energy, total_xs, elastic_xs, absorbtion_xs, fission_xs, nu_fission_xs)
    end

    # Close the file
    close(file)
    println("File closed successfully.")

    return immutableSimulationData(num_nucs, concs, mats, unionized_energy_array, index_grid, nuclide_grid, length_num_nucs, length_concs, length_mats, length_unionized_energy_array, length_index_grid, length_nuclide_grid, max_num_nucs, p_energy_samples, length_p_energy_samples, mat_samples, length_mat_samples)

end



function load_simulation_data()
    fname = "../openmp-offload/SimulationData_data.txt"
    println("Attempting to open file $fname...")
    println("The size of SimulationData is: ", sizeof(SimulationData))

    # Initialize SD
    SD = SimulationData(
        Vector{Int32}(undef, 0),
        Vector{Float64}(undef, 0),
        Vector{Int32}(undef, 0),
        Vector{Float64}(undef, 0),
        Vector{Int32}(undef, 0),
        Vector{NuclideGridPoint}(undef, 0),
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        Vector{Float64}(undef, 0),
        0,
        Vector{Int32}(undef, 0),
        0
    )

    # Attempt to open the file in read mode
    file = open(fname, "r")

    

    

    # Check if the file was successfully opened
    if file === nothing
        error("Failed to open file: $fname")
    end

    SD.length_num_nucs = 12
    SD.length_concs = 3852
    SD.length_mats = 3852
    SD.length_unionized_energy_array = 4012565
    SD.length_index_grid = 1424460575
    SD.length_nuclide_grid = 4012565
    SD.length_p_energy_samples = 0
    SD.length_mat_samples = 0
    SD.max_num_nucs = 321

    # Read the num_nucs value amount of integers from the file
    SD.num_nucs = Vector{Int32}(undef, SD.length_num_nucs)    

    # Read the num_nucs value amount of integers from the file
    for i in 1:SD.length_num_nucs
        SD.num_nucs[i] = parse(Int32, readline(file))
    end

    

    # Read the concs value amount of floats from the file
    SD.concs = Vector{Float64}(undef, SD.length_concs)

    # Read the concs value amount of floats from the file
    for i in 1:SD.length_concs
        SD.concs[i] = parse(Float64, readline(file))
    end

    # Read the mats value amount of integers from the file
    SD.mats = Vector{Int32}(undef, SD.length_mats)

    # Read the mats value amount of integers from the file
    for i in 1:SD.length_mats
        SD.mats[i] = parse(Int32, readline(file))
    end

    # Read the unionized_energy_array value amount of floats from the file
    SD.unionized_energy_array = Vector{Float64}(undef, SD.length_unionized_energy_array)

    # Read the unionized_energy_array value amount of floats from the file
    for i in 1:SD.length_unionized_energy_array
        SD.unionized_energy_array[i] = parse(Float64, readline(file))
    end

    # Read the index_grid value amount of integers from the file
    SD.index_grid = Vector{Int32}(undef, SD.length_index_grid)

    # Read the index_grid value amount of integers from the file
    for i in 1:SD.length_index_grid
        SD.index_grid[i] = parse(Int32, readline(file))
    end

    println("Got here")
    # Read the nuclide_grid value amount of NuclideGridPoint structs from the file
    SD.nuclide_grid = Vector{NuclideGridPoint}(undef, SD.length_nuclide_grid)

    # Initialize all the nuclide grids inside the SD
    for i in 1:SD.length_nuclide_grid
        SD.nuclide_grid[i] = NuclideGridPoint(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    end

    for i in 1:SD.length_nuclide_grid
        line = readline(file)
        SD.nuclide_grid[i].energy = parse(Float64, line)
        line = readline(file)
        SD.nuclide_grid[i].total_xs = parse(Float64, line)
        line = readline(file)
        SD.nuclide_grid[i].elastic_xs = parse(Float64, line)
        line = readline(file)
        SD.nuclide_grid[i].absorbtion_xs = parse(Float64, line)
        line = readline(file)
        SD.nuclide_grid[i].fission_xs = parse(Float64, line)
        line = readline(file)
        SD.nuclide_grid[i].nu_fission_xs = parse(Float64, line)
    end
   

    # Close the file
    close(file)
    println("File closed successfully.")
    
end

# Call the function to test it
# load_simulation_data()

# SD = load_sim_data()





end # module julia
