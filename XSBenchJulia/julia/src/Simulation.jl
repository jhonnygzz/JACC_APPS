module Simulation
include("LoadData.jl")
using .LoadData
import CUDA
using .CUDA
# using .julia
export run_event_based_simulation, hello
import JACC
# Grid types
const UNIONIZED = 0
const NUCLIDE = 1
const HASH = 2
const STARTING_SEED = 1070

function run_event_based_simulation(in:: LoadData.Input, SD:: LoadData.immutableSimulationData)

    println("Running baseline event-based simulation...")
    # Create JACC Arrays of all array data inside SD
    num_nucs_d = JACC.Array(SD.num_nucs)
    concs_d = JACC.Array(SD.concs)
    unionized_energy_array_d = JACC.Array(SD.unionized_energy_array)
    index_grid_d = JACC.Array(SD.index_grid)
    nuclide_grid_d = JACC.Array(SD.nuclide_grid)
    mats_d = JACC.Array(SD.mats)
    mat_samples_d = JACC.Array(SD.mat_samples)
    p_energy_samples_d = JACC.Array(SD.p_energy_samples)

    length_num_nucs = length(SD.num_nucs)
    length_concs = length(SD.concs)
    length_unionized_energy_array = length(SD.unionized_energy_array)
    length_index_grid = length(SD.index_grid)
    length_nuclide_grid = length(SD.nuclide_grid)
    length_mats = length(SD.mats)
    length_mat_samples = length(SD.mat_samples)
    length_p_energy_samples = length(SD.p_energy_samples)
    max_num_nucs = SD.max_num_nucs

    p_energy = 0.626011
    mat = 5
    n_iso = 355
    in_grid = 11303
    in_gridtype = 0
    in_hashbins = 10000
    macro_xs_vector = JACC.Array(zeros(Float64, 355))

    println("Lookups: ", in.lookups)
    # macro_xs_vector = CUDA.fill(Float64, 5)


    JACC.parallel_for(in.lookups, xs_lookup_kernel_baseline, in, num_nucs_d, concs_d,
     unionized_energy_array_d, index_grid_d, nuclide_grid_d, mats_d,
      mat_samples_d, p_energy_samples_d, length_num_nucs, length_concs,
      length_unionized_energy_array, length_index_grid, length_nuclide_grid,
      length_mats, length_mat_samples, length_p_energy_samples, max_num_nucs, p_energy, mat, n_iso, in_grid, in_gridtype, in_hashbins, macro_xs_vector)
    

    println("Simulation Complete.")
    println("LOL")
end





function xs_lookup_kernel_baseline(lookups, in, num_nucs_d, concs_d, unionized_energy_array_d,
     index_grid_d, nuclide_grid_d, mats_d, mat_samples_d, p_energy_samples_d, length_num_nucs,
     length_concs, length_unionized_energy_array, length_index_grid, length_nuclide_grid,
     length_mats, length_mat_samples, length_p_energy_samples, max_num_nucs, p_energy, mat,
      n_iso, in_grid, in_gridtype, in_hashbins, macro_xs_vector)
    # print("hi")
    # # Set the initial seed value
    # seed::UInt64 = STARTING_SEED

    # # Forward seed to lookup index (we need 2 samples per lookup)
	# seed = fast_forward_LCG(seed, 2*i); # TODO: Make fast_forward_LCG function

    # # Randomly pick an energy and material for the particle
    # p_energy = LCG_random_double(seed) # TODO: Make LCG_random_double function
    # mat = pick_mat(seed) # TODO: Make pick_mat function

    # a = 1 + 1

    # # Perform macroscopic Cross Section Lookup
    # calculate_macro_xs(p_energy, mat, length_num_nucs, length_unionized_energy_array, 
    # num_nucs_d, concs_d, unionized_energy_array_d, index_grid_d, 
    # nuclide_grid_d, mats_d, macro_xs_vector, UNIONIZED, 
    # length_index_grid, max_num_nucs)
    # add()
    # macro__xs()
    add()


end

function add()
    a = 1 + 1
end

function grid_search(n::Int64, quarry::Float64, A)::Int64
    lowerLimit = 1
    upperLimit = n
    length = upperLimit - lowerLimit

    while length > 1
        examinationPoint = lowerLimit + div(length, 2)
        
        if A[examinationPoint] > quarry
            upperLimit = examinationPoint
        else
            lowerLimit = examinationPoint
        end
        
        length = upperLimit - lowerLimit
    end
    
    return lowerLimit
end

function test(n::Int64, quarry::Float64, A::Vector{Float64})::Int64
    a = 1 + 1
end

function simulation_ex()
    println("Simulation Ex running Kernel")

    macro_xs_vector = JACC.Array(zeros(Float64, 355))
    grid_type = 0
    n_isotopes = 355
    n_gridpoints = 11303
    p_energy = 0.626011
    hash_bins = 10000
    mat = 5
    max_num_nucs = 5


    unionized_energy_array = JACC.Array(zeros(Float64, 11303)) # e_grid
    num_nucs = JACC.Array(zeros(Int64, 5))
    xs_vector = JACC.Array(zeros(Float64, 5))
    mats = JACC.Array(zeros(Int64, 5))
    concs = JACC.Array(zeros(Float64, 5))
    index_grid = JACC.Array(zeros(Int64, 355)) # index_data
    nu_grids = Vector{immutableNuclideGridPoint}()
    for i in 1:355
        push!(nu_grids, immutableNuclideGridPoint(0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    end
    nuclide_grid = JACC.Array(nu_grids)


    JACC.parallel_for(10, kernel, macro_xs_vector, grid_type, n_isotopes, n_gridpoints, p_energy, unionized_energy_array, hash_bins,
    mat, num_nucs, xs_vector, mats, max_num_nucs, concs, index_grid, nuclide_grid)
    println("Simulation Ex Complete.")
end

function kernel(lookups, macro_xs_vector, grid_type, n_isotopes, n_gridpoints, p_energy, unionized_energy_array, hash_bins,
    mat, num_nucs, xs_vector, mats, max_num_nucs, concs, index_grid, nuclide_grid)
    
    macro_xs(macro_xs_vector, grid_type, n_isotopes, n_gridpoints, p_energy, unionized_energy_array, hash_bins,
    mat, num_nucs, xs_vector, mats, max_num_nucs, concs, index_grid, nuclide_grid)
end

function macro_xs(macro_xs_vector, grid_type, n_isotopes, n_gridpoints, p_energy, egrid, hash_bins,
    mat, num_nucs, xs_vector, mats, max_num_nucs, concs, index_data, nuclide_grids)
    p_nuc = 0 # the nuclide we are looking up
    idx = -1
    conc = 0.0 # the concentration of the nuclide in the material

    # # cleans out macro_xs_vector
    for k in 1:5
        macro_xs_vector[k] = 0.0
    end

    if grid_type == UNIONIZED
        a = 1 + 1
        # idx = grid_search(n_isotopes * n_gridpoints, p_energy, egrid) TODO: Fix
    elseif grid_type == HASH
        du = 1.0 / hash_bins
        idx = Int(p_energy / du)
        end

    for j in 1:num_nucs[mat]
        # a = 1 + 1
        xs_vector = []
        p_nuc = mats[(mat-1)*max_num_nucs + (j-1) + 1]
        conc = concs[(mat-1)*max_num_nucs + (j-1) + 1]
        calculate_micro_xs(p_energy, p_nuc, n_isotopes, n_gridpoints, egrid, index_data,
        nuclide_grids, idx, xs_vector, grid_type, hash_bins)

    end
    



end

function calculate_macro_xs(p_energy::Float64, mat::Int, n_isotopes::Int64,
    n_gridpoints::Int64, num_nucs::Vector{Int64},
    concs::Vector{Float64}, egrid::Vector{Float64},
    index_data::Vector{Int}, nuclide_grids::Vector{NuclideGridPoint},
    mats::Vector{Int}, macro_xs_vector::Vector{Float64},
    grid_type::Int, hash_bins::Int, max_num_nucs::Int)

p_nuc = 0 # the nuclide we are looking up
idx = -1
conc = 0.0 # the concentration of the nuclide in the material

# cleans out macro_xs_vector
for k in 1:5
macro_xs_vector[k] = 0.0
end

# If we are using the unionized energy grid (UEG), we only
# need to perform 1 binary search per macroscopic lookup.
# If we are using the nuclide grid search, it will have to be
# done inside of the "calculate_micro_xs" function for each different
# nuclide in the material.
if grid_type == UNIONIZED
idx = grid_search(n_isotopes * n_gridpoints, p_energy, egrid)
elseif grid_type == HASH
du = 1.0 / hash_bins
idx = Int(p_energy / du)
end

# Once we find the pointer array on the UEG, we can pull the data
# from the respective nuclide grids, as well as the nuclide
# concentration data for the material
# Each nuclide from the material needs to have its micro-XS array
# looked up & interpolated (via calculate_micro_xs). Then, the
# micro XS is multiplied by the concentration of that nuclide
# in the material, and added to the total macro XS array.
# (Independent -- though if parallelizing, must use atomic operations
#  or otherwise control access to the xs_vector and macro_xs_vector to
#  avoid simultaneous writing to the same data structure)
for j in 1:num_nucs[mat]
xs_vector = zeros(Float64, 5)
p_nuc = mats[(mat-1)*max_num_nucs + (j-1) + 1]
conc = concs[(mat-1)*max_num_nucs + (j-1) + 1]
calculate_micro_xs(p_energy, p_nuc, n_isotopes, n_gridpoints, egrid, index_data,
   nuclide_grids, idx, xs_vector, grid_type, hash_bins)
for k in 1:5
macro_xs_vector[k] += xs_vector[k] * conc
end
end
end


function calculate_micro_xs(p_energy::Float64, nuc::Int, n_isotopes::Int64,
    n_gridpoints::Int64, egrid,
    index_data, nuclide_grids,
    idx::Int64, xs_vector, grid_type::Int, hash_bins::Int)

    # Variables
    f = 0.0
    low = immutableNuclideGridPoint(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    high = immutableNuclideGridPoint(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

        # If using only the nuclide grid, we must perform a binary search
        # to find the energy location in this particular nuclide's grid.
        if grid_type == NUCLIDE
            # Perform binary search on the Nuclide Grid to find the index
            idx = grid_search_nuclide(n_gridpoints, p_energy, view(nuclide_grids, (nuc*n_gridpoints + 1):(nuc*n_gridpoints + n_gridpoints)), 0, n_gridpoints-1)
            # pull ptr from nuclide grid and check to ensure that
            # we're not reading off the end of the nuclide's grid
            other = idx
            if idx == n_gridpoints - 1
                # low = nuclide_grids[nuc*n_gridpoints + idx - 1 + 1]
            else
                low = nuclide_grids[1]
                # low = nuclide_grids[nuc*n_gridpoints + idx + 1] # TODO LINE THAT BREAKS EVERYTHING
            end
        elseif grid_type == UNIONIZED
        a = 1 + 1
        index = idx * n_isotopes + nuc + 1
        data = index_data[index] 
        #     # Unionized Energy Grid - we already know the index, no binary search needed.
        #     # pull ptr from energy grid and check to ensure that
        #     # we're not reading off the end of the nuclide's grid
        #     index = idx * n_isotopes + nuc + 1
        #     # data = index_data[index]
        #     data = nuclide_grids[1]

        #     if 1 == n_gridpoints - 1
        #         # low = nuclide_grids[nuc*n_gridpoints + index_data[idx * n_isotopes + nuc + 1] - 1 + 1]
        #     else
        #         # low = nuclide_grids[nuc*n_gridpoints + index_data[idx * n_isotopes + nuc + 1] + 1]
        #     end
        end

end


function calculate_micro_xs_old(p_energy::Float64, nuc::Int, n_isotopes::Int64,
                            n_gridpoints::Int64, egrid,
                            index_data, nuclide_grids,
                            idx, xs_vector, grid_type::Int, hash_bins::Int)

    # Variables
    f = 0.0
    low = immutableNuclideGridPoint(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    high = immutableNuclideGridPoint(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # If using only the nuclide grid, we must perform a binary search
    # to find the energy location in this particular nuclide's grid.
    if grid_type == NUCLIDE
        # Perform binary search on the Nuclide Grid to find the index
        idx = grid_search_nuclide(n_gridpoints, p_energy, view(nuclide_grids, (nuc*n_gridpoints + 1):(nuc*n_gridpoints + n_gridpoints)), 0, n_gridpoints-1)
        # pull ptr from nuclide grid and check to ensure that
        # we're not reading off the end of the nuclide's grid
        if idx == n_gridpoints - 1
            low = nuclide_grids[nuc*n_gridpoints + idx - 1 + 1]
        else
            low = nuclide_grids[nuc*n_gridpoints + idx + 1]
        end

        

    elseif grid_type == UNIONIZED
        # Unionized Energy Grid - we already know the index, no binary search needed.
        # pull ptr from energy grid and check to ensure that
        # we're not reading off the end of the nuclide's grid
        if index_data[idx * n_isotopes + nuc + 1] == n_gridpoints - 1
            low = nuclide_grids[nuc*n_gridpoints + index_data[idx * n_isotopes + nuc + 1] - 1 + 1]
        else
            low = nuclide_grids[nuc*n_gridpoints + index_data[idx * n_isotopes + nuc + 1] + 1]
        end
    else
        # Hash grid
        # load lower bounding index
        u_low = index_data[idx * n_isotopes + nuc + 1]

        # Determine higher bounding index
        u_high = if idx == hash_bins - 1
            n_gridpoints - 1
        else
            index_data[(idx+1)*n_isotopes + nuc + 1] + 1
        end

        # Check edge cases to make sure energy is actually between these
        # Then, if things look good, search for gridpoint in the nuclide grid
        # within the lower and higher limits we've calculated.
        e_low  = nuclide_grids[nuc*n_gridpoints + u_low + 1].energy
        e_high = nuclide_grids[nuc*n_gridpoints + u_high + 1].energy
        lower = if p_energy <= e_low
            0
        elseif p_energy >= e_high
            n_gridpoints - 1
        else
            grid_search_nuclide(n_gridpoints, p_energy, nuclide_grids[nuc*n_gridpoints+1:nuc*n_gridpoints+n_gridpoints], u_low, u_high)
        end

        if lower == n_gridpoints - 1
            low = nuclide_grids[nuc*n_gridpoints + lower - 1 + 1]
        else
            low = nuclide_grids[nuc*n_gridpoints + lower + 1]
        end
    end

    high = nuclide_grids[nuc*n_gridpoints + (findfirst(==(low), nuclide_grids) + 1)]

    # calculate the re-useable interpolation factor
    f = (high.energy - p_energy) / (high.energy - low.energy)

    # Total XS
    xs_vector[1] = high.total_xs - f * (high.total_xs - low.total_xs)

    # Elastic XS
    xs_vector[2] = high.elastic_xs - f * (high.elastic_xs - low.elastic_xs)

    # Absorption XS
    xs_vector[3] = high.absorbtion_xs - f * (high.absorbtion_xs - low.absorbtion_xs)

    # Fission XS
    xs_vector[4] = high.fission_xs - f * (high.fission_xs - low.fission_xs)

    # Nu Fission XS
    xs_vector[5] = high.nu_fission_xs - f * (high.nu_fission_xs - low.nu_fission_xs)
end






function grid_search_nuclide(n::Int64, quarry::Float64, A, low::Int64, high::Int64)
    lowerLimit = low
    upperLimit = high
    length = upperLimit - lowerLimit

    while length > 1
        examinationPoint = lowerLimit + (length ÷ 2)
        
        if A[examinationPoint + 1].energy > quarry
            upperLimit = examinationPoint
        else
            lowerLimit = examinationPoint
        end
        
        length = upperLimit - lowerLimit
    end
    
    return lowerLimit
end





end # module

# function fast_forward_LCG(seed::UInt64, n::UInt64)::UInt64
#     # LCG parameters
#     const m = UInt64(9223372036854775808) # 2^63
#     a = UInt64(2806196910506780709)
#     c = UInt64(1)

#     n = n % m

#     a_new = UInt64(1)
#     c_new = UInt64(0)

#     while n > 0
#         if n & 1 != 0
#             a_new *= a
#             c_new = c_new * a + c
#         end
#         c *= (a + 1)
#         a *= a

#         n >>= 1
#     end

#     return (a_new * seed + c_new) % m
# end

