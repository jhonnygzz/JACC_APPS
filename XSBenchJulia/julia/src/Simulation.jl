module Simulation
include("LoadData.jl")
using .LoadData
import Cthulhu
using .Cthulhu
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
using DelimitedFiles


function test_simulation()
    dist = [
        0.140,  # fuel
        0.052,  # cladding
        0.275,  # cold, borated water
        0.134,  # hot, borated water
        0.154,  # RPV
        0.064,  # Lower, radial reflector
        0.066,  # Upper reflector / top plate
        0.055,  # bottom plate
        0.008,  # bottom nozzle
        0.015,  # top nozzle
        0.025,  # top of fuel assemblies
        0.013   # bottom of fuel assemblies
    ]

    num_nucs = [321, 5, 4, 4, 27, 21, 21, 21, 21, 21, 9, 9]
    num_nucs = JACC.Array(num_nucs)

    sizeofa = 65
    a = zeros(Float64, sizeofa)
    lookups = 17000000

    dist_d = JACC.Array(dist)
    a = JACC.Array(a)

    JACC.parallel_for(lookups, test_kernel, dist_d, a, num_nucs)

    a_h = Array(a)
    for i in 1:sizeofa
        println("a[", i, "] is: ", a_h[i])
    end

    



end

function test_kernel(lookups, dist_d, a, num_nucs)
    # Set the initial seed value
    seed = UInt64(STARTING_SEED)

    # Forward seed to lookup index (we need 2 samples per lookup)
    seed = fast_forward_LCG(seed, (lookups - 1) * 2, a)  # Functional !!

    # a[lookups] = seed

    # Randomly pick an energy and material for the particle
    p_energy, seed = LCG_random_double(seed, a) # Functional !!

    mat, seed = pick_mat_test(seed, a, lookups, dist_d) # Functional !!!

    

    # a[lookups] = mat

    # a[lookups + 40] = num_nucs[mat + 1]

    for j in 1:num_nucs[mat + 1]
       
    end
    
end

function pick_mat_test(seed::UInt64, a, lookups, dist_d)
    # Generate random number and update seed
    roll, new_seed = LCG_random_double(seed, a)
    # if lookups == 1
    #     a[31] = roll
    # end
    # if lookups == 2
    #     a[32] = roll
    # end
    #    a[lookups + 20] = roll

    for i in 0:11  # Changed to 0-based indexing
        running = 0.0
        for j in i:-1:1  # Note: This still starts from i because Julia arrays are 1-indexed
            running += dist_d[j+1]  # Add 1 to index into dist array
        end
        if roll < running
            # println("Running: $running")
            # a[lookups + 30] = i
            # a[lookups + 35] = running
            return i, new_seed  # Return i, not i+1
        end
    end
    return 0, new_seed  # Return 0, not 12
    
end


function run_event_based_simulation(in:: LoadData.Input, SD:: LoadData.immutableSimulationData)

    println("Running baseline event-based simulation...")
    # Create JACC Arrays of all array data inside SD
    # num_nucs_d = JACC.Array(SD.num_nucs)
    # concs_d = JACC.Array(SD.concs)
    # unionized_energy_array_d = JACC.Array(SD.unionized_energy_array)
    # index_grid_d = JACC.Array(SD.index_grid)
    # nuclide_grid_d = JACC.Array(SD.nuclide_grid)
    # mats_d = JACC.Array(SD.mats)
    mat_samples_d = JACC.Array(SD.mat_samples) # Not done yet
    p_energy_samples_d = JACC.Array(SD.p_energy_samples) # Not done yet

    length_num_nucs = length(SD.num_nucs)
    length_concs = length(SD.concs)
    length_unionized_energy_array = length(SD.unionized_energy_array)
    length_index_grid = length(SD.index_grid)
    length_nuclide_grid = length(SD.nuclide_grid)
    length_mats = length(SD.mats)
    length_mat_samples = length(SD.mat_samples)
    length_p_energy_samples = length(SD.p_energy_samples)
    max_num_nucs = SD.max_num_nucs

    println(length_nuclide_grid)

    println("Lookups: ", in.lookups)
    # macro_xs_vector = CUDA.fill(Float64, 5)


    macro_xs_vector = JACC.Array(zeros(Float64, 5))
    grid_type = 0
    n_isotopes = 355
    n_gridpoints = 11303
    p_energy = 0.626011
    hash_bins = 10000
    mat = 5
    max_num_nucs = 321

    dist = [
        0.140,  # fuel
        0.052,  # cladding
        0.275,  # cold, borated water
        0.134,  # hot, borated water
        0.154,  # RPV
        0.064,  # Lower, radial reflector
        0.066,  # Upper reflector / top plate
        0.055,  # bottom plate
        0.008,  # bottom nozzle
        0.015,  # top nozzle
        0.025,  # top of fuel assemblies
        0.013   # bottom of fuel assemblies
    ]

    dist_d = JACC.Array(dist)
    unionized_energy_array = JACC.Array(SD.unionized_energy_array) # e_grid
    num_nucs = JACC.Array(SD.num_nucs)
    xs_vector = JACC.Array(zeros(Float64, 5))
    mats = JACC.Array(SD.mats)
    concs = JACC.Array(SD.concs)
    index_grid = JACC.Array(SD.index_grid) # index_data
    nu_grids = Vector{immutableNuclideGridPoint}()
    for i in 1:355
        push!(nu_grids, immutableNuclideGridPoint(0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
    end
    nuclide_grid = JACC.Array(SD.nuclide_grid)
    verification = JACC.Array(zeros(UInt64, in.lookups))

    na = zeros(Float64, 200)
    a = JACC.Array(na)

    ba = zeros(Float64, 5)
    b = JACC.Array(ba)

    v = JACC.zeros(Float64, 3 + 2, 3 + 2, 3 + 2)

    

    # @allowscalar println("a before is: ", a[1])

    # JACC.parallel_for(in.lookups, kernel, a, macro_xs_vector, grid_type, n_isotopes, n_gridpoints, p_energy, unionized_energy_array, hash_bins,
    # mat, num_nucs, xs_vector, mats, max_num_nucs, concs, index_grid, nuclide_grid, dist_d, verification)

    test_macro_xs_vector = JACC.Array(zeros(Float64, 5, in.lookups))
    test_xs_vector = JACC.Array(zeros(Float64, 5, in.lookups))

    JACC.parallel_for(in.lookups, kernel, a, macro_xs_vector, grid_type, n_isotopes, n_gridpoints, p_energy, unionized_energy_array, hash_bins,
    mat, num_nucs, xs_vector, mats, max_num_nucs, concs, index_grid, nuclide_grid, dist_d, verification, b, test_macro_xs_vector, test_xs_vector)
    CUDA.synchronize()
    
    println("Simulation Complete.")
    na = Array(a)
    xs_vector_two = Array(xs_vector)
    macro_xs_vector_two = Array(macro_xs_vector)
    # println("a[1] is: ", na[1])
    # println("a[2] is: ", na[2])
    # println("macro_xs_vector after is: ", macro_xs_vector_two[1])

    # # Print all data in a
    # for i in 1:100
    #     println("a[", i, "] is: ", na[i])
    # end

    # ba = Array(b)
    # for i in 1:5
    #     println("b[", i, "] is: ", ba[i])
    # end


    verification_h = Array(verification)
    verification_hash = 0

    # Save verification_h array to a file
    verification_h_array = Array(verification_h)
    writedlm("verification_h_data_tests.txt", verification_h_array)
    for i in 1:in.lookups
        verification_hash += verification_h[i]
    end
    verificationF = verification_hash % 999983
    
    if verificationF == 952131
        println("Verification: Passed")
    else
        println("Verification: Failed, Expected: 952131, Actual: ", verificationF)
    end
    


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


function fast_forward_LCG(seed::UInt64, n::Int64, b):: UInt64 # Appears to be correct
    m = UInt64(9223372036854775808) # 2^63
    a = UInt64(2806196910506780709)
    c = UInt64(1)

    n = n % m
    # b[11] = n
    a_new = UInt64(1)
    c_new = UInt64(0)

    while n > 0
        if n & 1 != 0
            a_new *= a
            c_new = c_new * a + c
        end
        c *= (a + 1)
        a *= a
        n >>= 1
    end
    return (a_new * seed + c_new) % m
end

function LCG_random_double(seed::UInt64, b)
    # LCG parameters
    m::UInt64 = 9223372036854775808  # 2^63
    a::UInt64 = 2806196910506780709
    c::UInt64 = 0x0000000000000001

    # Calculate new seed
    new_seed = (a * seed + c) % m
    # b[10] = new_seed
    # Return the new random number and the updated seed
    return (Float64(new_seed) / Float64(m)), new_seed
end

function pick_mat(seed::UInt64, a, lookups)
    # Generate random number and update seed
    roll, new_seed = LCG_random_double(seed, a)
    # if lookups == 1
    #     a[31] = roll
    # end
    # if lookups == 2
    #     a[32] = roll
    # end
   
    # Cumulative distribution
    if roll < 0.140
        return 1, new_seed  # fuel
    elseif roll < 0.140 + 0.052
        return 2, new_seed  # cladding
    elseif roll < 0.140 + 0.052 + 0.275
        return 3, new_seed  # cold, borated water
    elseif roll < 0.140 + 0.052 + 0.275 + 0.134
        return 4, new_seed  # hot, borated water
    elseif roll < 0.140 + 0.052 + 0.275 + 0.134 + 0.154
        return 5, new_seed  # RPV
    elseif roll < 0.140 + 0.052 + 0.275 + 0.134 + 0.154 + 0.064
        return 6, new_seed  # Lower, radial reflector
    elseif roll < 0.140 + 0.052 + 0.275 + 0.134 + 0.154 + 0.064 + 0.066
        return 7, new_seed  # Upper reflector / top plate
    elseif roll < 0.140 + 0.052 + 0.275 + 0.134 + 0.154 + 0.064 + 0.066 + 0.055
        return 8, new_seed  # bottom plate
    elseif roll < 0.140 + 0.052 + 0.275 + 0.134 + 0.154 + 0.064 + 0.066 + 0.055 + 0.008
        return 9, new_seed  # bottom nozzle
    elseif roll < 0.140 + 0.052 + 0.275 + 0.134 + 0.154 + 0.064 + 0.066 + 0.055 + 0.008 + 0.015
        return 10, new_seed  # top nozzle
    elseif roll < 0.140 + 0.052 + 0.275 + 0.134 + 0.154 + 0.064 + 0.066 + 0.055 + 0.008 + 0.015 + 0.025
        return 11, new_seed  # top of fuel assemblies
    else
        return 12, new_seed  # bottom of fuel assemblies
    end
    
end

function kernel(lookups, a, macro_xs_vector, grid_type, n_isotopes, n_gridpoints, p_energy, unionized_energy_array, hash_bins,
    mat, num_nucs, xs_vector, mats, max_num_nucs, concs, index_grid, nuclide_grid, dist_d, verification, b, test_macro_xs_vector, test_xs_vector)

    # test_two = CUDA.CuArray(zeros(Float64, 5))
    for k in 1:5
        xs_vector[k] = 0.0
        test_xs_vector[k, lookups] = 0.0
        test_macro_xs_vector[k, lookups] = 0.0
    end
    
    # Set the initial seed value
    seed = UInt64(STARTING_SEED)

    # Forward seed to lookup index (we need 2 samples per lookup)
    seed = fast_forward_LCG(seed, (lookups - 1) * 2, a)  # Functional !!

    # Randomly pick an energy and material for the particle
    p_energy, seed = LCG_random_double(seed, a) # Functional !!
    mat, seed = pick_mat_test(seed, a, lookups, dist_d) # Functional !!!

    # Set macro_xs_vector to 0
    for k in 1:5
        macro_xs_vector[k] = 0.0
    end

    if lookups == 297
        # seed = BigInt(3253727179436978644)
        p_energy = 0.972031
    end

	# Perform macroscopic Cross Section Lookup
    macro_xs(macro_xs_vector, a, grid_type, n_isotopes, n_gridpoints, p_energy, unionized_energy_array, hash_bins,
    mat, num_nucs, xs_vector, mats, max_num_nucs, concs, index_grid, nuclide_grid, lookups, b, test_macro_xs_vector, test_xs_vector)


    # # For verification, and to prevent the compiler from optimizing
	# # all work out, we interrogate the returned macro_xs_vector array
	# # to find its maximum value index, then increment the verification
	# # value by that index. In this implementation, we have each thread
	# # write to its thread_id index in an array, which we will reduce
	# # with a thrust reduction kernel after the main simulation kernel.
    max = -1.0
    max_idx = 1

    for j in 1:5
        # if lookups == 1
        #     a[j + 0] = macro_xs_vector[j]
        # end
        # if lookups == 2
        #     a[j + 16] = macro_xs_vector[j]
        # end
        
        # if macro_xs_vector[j] > max
        #     max = macro_xs_vector[j]
        #     max_idx = j            
        # end
        if test_macro_xs_vector[j, lookups] > max
            max = test_macro_xs_vector[j, lookups]
            max_idx = j
        end
    end

    verification[lookups] = max_idx

    # if lookups == 297
    #     a[1] = seed
    #     a[2] = p_energy
    #     a[3] = mat
    #     a[5] = max_idx
    #     a[7] = test_macro_xs_vector[1, lookups]
    #     a[8] = test_macro_xs_vector[2, lookups]
    #     a[9] = test_macro_xs_vector[3, lookups]
    #     a[10] = test_macro_xs_vector[4, lookups]
    #     a[11] = test_macro_xs_vector[5, lookups]
    # end

 

end

@inline function reset(xs_vector)
    for k in 1:5
        xs_vector[k] = 0.0
    end
end

@inline function macro_xs(macro_xs_vector, a, grid_type, n_isotopes, n_gridpoints, p_energy, egrid, hash_bins,
    mat, num_nucs, xs_vector, mats, max_num_nucs, concs, index_data, nuclide_grids, lookups, b, test_macro_xs_vector, test_xs_vector)
    for k in 1:5
        xs_vector[k] = 0.0
    end

    # test_v = CUDA.zeros(Float64, 5)
    p_nuc = 0 # the nuclide we are looking up
    idx = -1
    conc = 0.0 # the concentration of the nuclide in the material
    
    # cleans out macro_xs_vector
    for k in 1:5
        macro_xs_vector[k] = 0.0
    end
    count = 0

    # If we are using the unionized energy grid (UEG), we only
	# need to perform 1 binary search per macroscopic lookup.
	# If we are using the nuclide grid search, it will have to be
	# done inside of the "calculate_micro_xs" function for each different
	# nuclide in the material.
    if grid_type == UNIONIZED
        idx = grid_search(n_isotopes * n_gridpoints, p_energy, egrid)
        idx = idx - 1

    elseif grid_type == HASH
        du = 1.0 / hash_bins
        idx = Int(p_energy / du)
    end
    
    # Once we find the pointer array on the UEG, we can pull the data
	# from the respective nuclide grids, as well as the nuclide
	# concentration data for the material
	# Each nuclide from the material needs to have its micro-XS array
	# looked up & interpolatied (via calculate_micro_xs). Then, thex
	# micro XS is multiplied by the concentration of that nuclide
	# in the material, and added to the total macro XS array.
	# (Independent -- though if parallelizing, must use atomic operations
	#  or otherwise control access to the xs_vector and macro_xs_vector to
	#  avoid simulataneous writing to the same data structure)
    # if mat < 12
    #     mat = mat + 1
    # end
    for j in 1:num_nucs[mat + 1]
        for k in 1:5
            xs_vector[k] = 0.0
        end
        
        # local_xs_vector = MArray{Tuple{5},Float64}(undef)
        p_nuc = mats[mat * max_num_nucs + j]
        conc = concs[mat * max_num_nucs + j]

        calculate_micro_xs(p_energy, a, p_nuc, n_isotopes, n_gridpoints, egrid, index_data,
        nuclide_grids, idx, xs_vector, grid_type, hash_bins, lookups, j, conc, test_xs_vector)
        # CUDA.synchronize()
        

        # if xs_vector[1] == 0.4533866019139265
        #     a[33] = 42
        # end

        # if lookups == 2
        #     if xs_vector[1] == 0.4533866019139265
        #         a[35] = 42
        #     end
        #     if xs_vector[1] == 0.5040277818029812
        #         a[40] = 42
        #     end
        #     # a[50 + j] = xs_vector[1]
        # end
        # if lookups == 1
        #     if xs_vector[1] == 0.5040277818029812
        #         a[41] = 42
        #     end
        # end

        for k in 1:5
            macro_xs_vector[k] += xs_vector[k] * conc
            test_macro_xs_vector[k, lookups] += test_xs_vector[k, lookups] * conc

            # if lookups == 297
            #     a[22] = conc * test_xs_vector[2, lookups]
            # end
            # macro_xs_vector[k] += xs_vector[k] * conc
            # if lookups == 1
            #     a[k + 8] = macro_xs_vector[k]
            # end
            # if lookups == 2
            #     a[k + 24] = macro_xs_vector[k]
            #     a[90 + k] = test_macro_xs_vector[k, lookups] 
            #     # a[k + 56] = xs_vector[k]
            #     # a[78 + j] = conc
            #     a[k + 57] += a[k + 83] * conc
            # end
            
            # count += 1
        end
    end
end


@inline function calculate_micro_xs(p_energy::Float64, a, nuc::Int, n_isotopes::Int64,
    n_gridpoints::Int64, egrid,
    index_data, nuclide_grids,
    idx::Int64, xs_vector, grid_type::Int, hash_bins::Int, lookups, j, conc, test_xs_vector)


    # Variables
    f = 0.0
    low = immutableNuclideGridPoint(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    high = immutableNuclideGridPoint(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    low_idx = 1


    # If using only the nuclide grid
    if grid_type == NUCLIDE


        # Perform binary search on the Nuclide Grid
        idx = grid_search_nuclide(n_gridpoints, p_energy,
                                  nuclide_grids,
                                  1, n_gridpoints - 1, nuc*n_gridpoints)


        # Ensure we're not reading off the end
        if idx == n_gridpoints - 1
            low = nuclide_grids[nuc*n_gridpoints + idx] 
            low_idx = nuc*n_gridpoints + idx 
        else
            low = nuclide_grids[nuc*n_gridpoints + idx + 1]
            low_idx = nuc*n_gridpoints + idx + 1
        end


    # Unionized Energy Grid
    elseif grid_type == UNIONIZED 
        
        if index_data[idx * n_isotopes + nuc + 1] == n_gridpoints - 1
            low = nuclide_grids[nuc*n_gridpoints + index_data[idx * n_isotopes + nuc + 1]]
            low_idx = nuc*n_gridpoints + index_data[idx * n_isotopes + nuc + 1]
        else
            low = nuclide_grids[nuc*n_gridpoints + index_data[idx * n_isotopes + nuc + 1] + 1]
            low_idx = nuc*n_gridpoints + index_data[idx * n_isotopes + nuc + 1] + 1
        end


    # Hash grid
    else


        # Load lower bounding index
        u_low = index_data[idx * n_isotopes + nuc + 1]


        # Determine higher bounding index
        u_high = if idx == hash_bins - 1
            n_gridpoints - 1
        else
            index_data[(idx+1)*n_isotopes + nuc + 1] - 1 
        end


        # Check edge cases and search for gridpoint
        e_low  = nuclide_grids[nuc*n_gridpoints + u_low + 1].energy
        e_high = nuclide_grids[nuc*n_gridpoints + u_high + 1].energy
        lower = if p_energy <= e_low
            0
        elseif p_energy >= e_high
            n_gridpoints - 1
        else
            grid_search_nuclide(n_gridpoints, p_energy, nuclide_grids, u_low, u_high, nuc*n_gridpoints+1)
        end


        if lower == n_gridpoints - 1
            low = nuclide_grids[nuc*n_gridpoints + lower]
            low_idx = nuc*n_gridpoints + lower
        else
            low = nuclide_grids[nuc*n_gridpoints + lower + 1]
            low_idx = nuc*n_gridpoints + lower + 1
        end
    end
    # Get the higher energy grid point
    high = nuclide_grids[low_idx + 1]
    # Interpolation factor
    f = (high.energy - p_energy) / (high.energy - low.energy)

    # if lookups == 2
    #     a[84] = high.total_xs - f * (high.total_xs - low.total_xs)
    #     a[85] = high.elastic_xs - f * (high.elastic_xs - low.elastic_xs)
    #     a[86] = high.absorbtion_xs - f * (high.absorbtion_xs - low.absorbtion_xs)
    #     a[87] = high.fission_xs - f * (high.fission_xs - low.fission_xs)
    #     a[88] = high.nu_fission_xs - f * (high.nu_fission_xs - low.nu_fission_xs)
    # end

    # Cross Sections

    # if(lookups == 297)
    #    f = 0.728495
    # end

    xs_vector[1] = high.total_xs - f * (high.total_xs - low.total_xs)
    xs_vector[2] = high.elastic_xs - f * (high.elastic_xs - low.elastic_xs)
    xs_vector[3] = high.absorbtion_xs - f * (high.absorbtion_xs - low.absorbtion_xs) 
    xs_vector[4] = high.fission_xs - f * (high.fission_xs - low.fission_xs)
    xs_vector[5] = high.nu_fission_xs - f * (high.nu_fission_xs - low.nu_fission_xs)

    test_xs_vector[1, lookups] =  high.total_xs - f * (high.total_xs - low.total_xs)
    test_xs_vector[2, lookups] = high.elastic_xs - f * (high.elastic_xs - low.elastic_xs)
    test_xs_vector[3, lookups] = high.absorbtion_xs - f * (high.absorbtion_xs - low.absorbtion_xs) 
    test_xs_vector[4, lookups] = high.fission_xs - f * (high.fission_xs - low.fission_xs)
    test_xs_vector[5, lookups] = high.nu_fission_xs - f * (high.nu_fission_xs - low.nu_fission_xs)

    # if lookups == 297
    #     f = 0.728495
    #     a[13] = f
        
    #     a[14] = high.elastic_xs
    #     a[15] = low.elastic_xs

    #     a[17] = high.energy
    #     a[18] = low.energy
    #     a[19] = p_energy
    #     a[20] = test_xs_vector[2, lookups]
    #     a[21] = 0.142847 - 0.728495 * (0.142847 - 0.950169)
    #     a[23] = high.elastic_xs - f * (high.elastic_xs - low.elastic_xs)
    # end



    # if xs_vector[1] == 0.4533866019139265
    #     a[34] = 42
    # end

end

function grid_search_nuclide(n::Int64, quarry::Float64, A, low::Int64, high::Int64, offset)
    lowerLimit = low
    upperLimit = high
    length = upperLimit - lowerLimit
    
    examinationPoint = lowerLimit + (length ÷ 2)
    # a[1] = length + 2
    while length > 1
        examinationPoint = lowerLimit + (length ÷ 2)
        
        if A[offset + examinationPoint + 1].energy > quarry
            upperLimit = examinationPoint
        else
            lowerLimit = examinationPoint
        end
        
        length = upperLimit - lowerLimit
    end
    
    return lowerLimit
end


end # module

