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
    verification = JACC.Array(zeros(UInt64, 17000001))

    na = zeros(Float64, 200)
    a = JACC.Array(na)

    ba = zeros(Float64, 5)
    b = JACC.Array(ba)


    
    # @allowscalar println("a before is: ", a[1])

    # JACC.parallel_for(in.lookups, kernel, a, macro_xs_vector, grid_type, n_isotopes, n_gridpoints, p_energy, unionized_energy_array, hash_bins,
    # mat, num_nucs, xs_vector, mats, max_num_nucs, concs, index_grid, nuclide_grid, dist_d, verification)

    JACC.parallel_for(in.lookups, kernel, a, macro_xs_vector, grid_type, n_isotopes, n_gridpoints, p_energy, unionized_energy_array, hash_bins,
    mat, num_nucs, xs_vector, mats, max_num_nucs, concs, index_grid, nuclide_grid, dist_d, verification, b)
    
    println("Simulation Complete.")
    na = Array(a)
    xs_vector_two = Array(xs_vector)
    macro_xs_vector_two = Array(macro_xs_vector)
    # println("a[1] is: ", na[1])
    # println("a[2] is: ", na[2])
    # println("macro_xs_vector after is: ", macro_xs_vector_two[1])

    # Print all data in a
    for i in 1:200
        println("a[", i, "] is: ", na[i])
    end


    verification_h = Array(verification)
    verification_hash = 0
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
    # LCG parameters
    # b[6] = b[6] + 1
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

function pick_mat(seed::UInt64, a)
    # Generate random number and update seed
    roll, new_seed = LCG_random_double(seed, a)
    # a[20] = roll
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

function add(seed::UInt64, n::Int64)::UInt64
    return 1 + 1
end

function kernel(lookups, a, macro_xs_vector, grid_type, n_isotopes, n_gridpoints, p_energy, unionized_energy_array, hash_bins,
    mat, num_nucs, xs_vector, mats, max_num_nucs, concs, index_grid, nuclide_grid, dist_d, verification, b)
    
    # Set the initial seed value
    seed = UInt64(STARTING_SEED)

    # Forward seed to lookup index (we need 2 samples per lookup)
    seed = fast_forward_LCG(seed, (lookups - 1) * 2, a)  # Functional !!

    a[6] = seed

    a[21] = lookups
    # a[6] = seed
    if lookups == 8
        a[6] = 8
    end
    # Randomly pick an energy and material for the particle
    p_energy, seed = LCG_random_double(seed, a) # Functional !!
    
    # mat, seed = pick_mat(seed, dist_d, a) 
    mat, seed = pick_mat(seed, a) # Functional !!!

    a[10] = p_energy
    a[11] = mat

    # a[7] = p_energy 
    # a[8] = mat 
    # Set macro_xs_vector to 0
    for k in 1:5
        macro_xs_vector[k] = 0.0
    end

    # macro_xs_vector = []
    # a[1] = mat
    macro_xs(macro_xs_vector, a, grid_type, n_isotopes, n_gridpoints, p_energy, unionized_energy_array, hash_bins,
    mat, num_nucs, xs_vector, mats, max_num_nucs, concs, index_grid, nuclide_grid, lookups, b)

    max = -1.0
    max_idx = 1

    for j in 1:5
        if lookups == 2
            a[j + 20] = macro_xs_vector[j]
        end
        if lookups == 1
            a[j + 30] = macro_xs_vector[j]
        end
        if macro_xs_vector[j] > max
            max = macro_xs_vector[j]
            max_idx = j
            
        end
    end


    verification[lookups] = max_idx

end

function macro_xs(macro_xs_vector, a, grid_type, n_isotopes, n_gridpoints, p_energy, egrid, hash_bins,
    mat, num_nucs, xs_vector, mats, max_num_nucs, concs, index_data, nuclide_grids, lookups, b)
    p_nuc = 0 # the nuclide we are looking up
    idx = -1
    conc = 0.0 # the concentration of the nuclide in the material
    
    # cleans out macro_xs_vector
    for k in 1:5
        macro_xs_vector[k] = 0.0
    end

    if grid_type == UNIONIZED
        idx = grid_search(n_isotopes * n_gridpoints, p_energy, egrid)
        idx = idx - 1

        a[12] = n_isotopes * n_gridpoints
        a[13] = p_energy
        a[14] = idx
    elseif grid_type == HASH
        du = 1.0 / hash_bins
        idx = Int(p_energy / du)
        a[13] = 55
    end
    if lookups == 2
        a[15] = num_nucs[mat + 1]
    end
   
    # if lookups == 2
    #     # a[52] = mat + 1
    #     # a[51] = num_nucs[mat + 1]
    # end
    
    for j in 1:num_nucs[mat + 1]
        
        for k in 1:5
            xs_vector[k] = 0.0
        end
        # p_nuc = mats[(mat-1)*max_num_nucs + (j-1) + 1]
        # conc = concs[(mat-1)*max_num_nucs + (j-1) + 1]
        p_nuc = mats[mat * max_num_nucs + j]
        
        conc = concs[mat * max_num_nucs + j]
        a[99] = (mat-1) # *max_num_nucs + (j-1) + 1
        a[18] = max_num_nucs
        # a[6] =  (j-1)
        a[7] = mat
        a[8] = j

        a[19] = p_nuc
        a[20] = conc
        # p_nuc = mats[mat * max_num_nucs + j]
        # conc = concs[mat * max_num_nucs + j]
        # if lookups == 2
        #     for k in 1:5
        #         xs_vector[k] = 0.0
        #         a[100] = 10000
        #     end
        # end

        if lookups == 2
            if xs_vector[1] == 0
                a[115] = 42
            end            
        end

        for k in 1:5
            b[k] = 0.0
        end


        calculate_micro_xs(p_energy, a, p_nuc, n_isotopes, n_gridpoints, egrid, index_data,
        nuclide_grids, idx, xs_vector, grid_type, hash_bins, lookups, j, conc, b)
        if lookups == 2
            if xs_vector[1] == 0.4533866019139265
                a[102] = 42
            end
            if xs_vector[1] == 0.5040277818029812
                a[110] = 42
            end
            a[105 + j] = xs_vector[1]
        end
        if lookups == 1
            if xs_vector[1] == 0.5040277818029812
                a[101] = 42
            end
        end
        for k in 1:5
            # macro_xs_vector[k] += xs_vector[k] * conc
            macro_xs_vector[k] += xs_vector[k] * conc
            a[194 + k] += b[k] * conc
            if lookups == 1
                a[169 + k] = a[194 + k]
                a[k + 38] = macro_xs_vector[k]
                a[k + 43] = xs_vector[k]
                a[49] = conc
            end
            if lookups == 2
                a[175 + k] = a[194 + k]
                a[k + 51] = macro_xs_vector[k]
                # a[88] += a[j + 73] * conc
                a[k + 56] = xs_vector[k]
                a[78 + j] = conc
                a[k + 88] += a[k + 83] * conc
                if j == 4 && k == 5
                    for b in 1:5
                        a[99] += a[b + 73] * a[78 + b]
                    end
                end
            end
        end
    end
end




function calculate_micro_xs(p_energy::Float64, a, nuc::Int, n_isotopes::Int64,
    n_gridpoints::Int64, egrid,
    index_data, nuclide_grids,
    idx::Int64, xs_vector, grid_type::Int, hash_bins::Int, lookups, j, conc, b)


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


        if lookups == 2
            a[26] = index_data[idx * n_isotopes + nuc + 1]
            a[27] = n_gridpoints - 1
        end
        
        if index_data[idx * n_isotopes + nuc + 1] == n_gridpoints - 1
            low = nuclide_grids[nuc*n_gridpoints + index_data[idx * n_isotopes + nuc + 1]]
            low_idx = nuc*n_gridpoints + index_data[idx * n_isotopes + nuc + 1]
        else
            low = nuclide_grids[nuc*n_gridpoints + index_data[idx * n_isotopes + nuc + 1] + 1]
            low_idx = nuc*n_gridpoints + index_data[idx * n_isotopes + nuc + 1] + 1


            if lookups == 1 && j == 4
                a[28] = low_idx
            end
            if lookups == 2 && j == 4
                a[29] = low_idx
            end
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

    if lookups == 2 && j == 4
        a[37] = f
    end
    a[9] = a[9] + 1
    
    a[16] = high.energy - p_energy
    a[17] = high.energy - low.energy


    if lookups == 1
        a[j + 63] = high.total_xs - f * (high.total_xs - low.total_xs)
        if j == 5
            a[95] += a[j + 63] * conc
            a[96] = conc
        end
    end
    if lookups == 2
        a[j + 73] = high.total_xs - f * (high.total_xs - low.total_xs)
        if j == 4
            a[97] += a[j + 73] * conc
            a[98] = conc
        end
    end

    if lookups == 2
        a[84] = high.total_xs - f * (high.total_xs - low.total_xs)
        a[85] = high.elastic_xs - f * (high.elastic_xs - low.elastic_xs)
        a[86] = high.absorbtion_xs - f * (high.absorbtion_xs - low.absorbtion_xs)
        a[87] = high.fission_xs - f * (high.fission_xs - low.fission_xs)
        a[88] = high.nu_fission_xs - f * (high.nu_fission_xs - low.nu_fission_xs)
    end



    # Cross Sections
    xs_vector[1] = high.total_xs - f * (high.total_xs - low.total_xs)
    xs_vector[2] = high.elastic_xs - f * (high.elastic_xs - low.elastic_xs)
    xs_vector[3] = high.absorbtion_xs - f * (high.absorbtion_xs - low.absorbtion_xs) 
    xs_vector[4] = high.fission_xs - f * (high.fission_xs - low.fission_xs)
    xs_vector[5] = high.nu_fission_xs - f * (high.nu_fission_xs - low.nu_fission_xs)

    b[1] = high.total_xs - f * (high.total_xs - low.total_xs)
    b[2] = high.elastic_xs - f * (high.elastic_xs - low.elastic_xs)
    b[3] = high.absorbtion_xs - f * (high.absorbtion_xs - low.absorbtion_xs) 
    b[4] = high.fission_xs - f * (high.fission_xs - low.fission_xs)
    b[5] = high.nu_fission_xs - f * (high.nu_fission_xs - low.nu_fission_xs)

    if xs_vector[1] == 0.4533866019139265
        a[105] = 42
    end


end

# function calculate_micro_xs(p_energy::Float64, a, nuc::Int, n_isotopes::Int64,
#     n_gridpoints::Int64, egrid,
#     index_data, nuclide_grids,
#     idx::Int64, xs_vector, grid_type::Int, hash_bins::Int, lookups, j, conc)

#     # # Variables
#     f = 0.0
#     low = immutableNuclideGridPoint(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
#     high = immutableNuclideGridPoint(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
#     low_idx = 1

#         # If using only the nuclide grid, we must perform a binary search
#         # to find the energy location in this particular nuclide's grid.
#         if grid_type == NUCLIDE# Originally NUCLIDE

#             # Perform binary search on the Nuclide Grid to find the index
#             idx = grid_search_nuclide(n_gridpoints, p_energy, 
#                           nuclide_grids, 
#                           1, n_gridpoints - 1, nuc*n_gridpoints)
#             # pull ptr from nuclide grid and check to ensure that
#             # we're not reading off the end of the nuclide's grid
#             if idx == n_gridpoints - 1
#                 low = nuclide_grids[nuc*n_gridpoints + idx - 1 + 1]
#                 low_idx = nuc*n_gridpoints + idx - 1 + 1
#             else
#                 low = nuclide_grids[nuc*n_gridpoints + idx + 1]
#                 low_idx = nuc*n_gridpoints + idx + 1
#             end
#         elseif grid_type == UNIONIZED # Unionized Energy Grid - we already know the index, no binary search needed.
#             # pull ptr from energy grid and check to ensure that
#             # we're not reading off the end of the nuclide's grid
#             if lookups == 2
#                 a[26] = index_data[idx * n_isotopes + nuc + 1]
#                 a[27] = n_gridpoints - 1
#             end
            
#             if index_data[idx * n_isotopes + nuc + 1] == n_gridpoints - 1
#                 low = nuclide_grids[nuc*n_gridpoints + index_data[idx * n_isotopes + nuc + 1] - 1 + 1]
#                 low_idx = nuc*n_gridpoints + index_data[idx * n_isotopes + nuc + 1] - 1 + 1
#             else
#                 low = nuclide_grids[nuc*n_gridpoints + index_data[idx * n_isotopes + nuc + 1] + 1]
#                 low_idx = nuc*n_gridpoints + index_data[idx * n_isotopes + nuc + 1] + 1
#                 if lookups == 1 && j == 4
#                     a[28] = low_idx
#                 end
#                 if lookups == 2 && j == 4
#                     a[29] = low_idx
#                 end

#             end
#         else
#             # Hash grid
#             # load lower bounding index
#             u_low = index_data[idx * n_isotopes + nuc + 1]

#             # Determine higher bounding index
#             u_high = if idx == hash_bins - 1
#                 n_gridpoints - 1
#             else
#                 index_data[(idx+1)*n_isotopes + nuc + 1] + 1
#             end

#             # Check edge cases to make sure energy is actually between these
#             # Then, if things look good, search for gridpoint in the nuclide grid
#             # within the lower and higher limits we've calculated.
#             e_low  = nuclide_grids[nuc*n_gridpoints + u_low + 1].energy
#             e_high = nuclide_grids[nuc*n_gridpoints + u_high + 1].energy
#             lower = if p_energy <= e_low
#                 0
#             elseif p_energy >= e_high
#                 n_gridpoints - 1
#             else
#                 grid_search_nuclide(n_gridpoints, p_energy, nuclide_grids, u_low, u_high, nuc*n_gridpoints+1)
#             end

#             if lower == n_gridpoints - 1
#                 low = nuclide_grids[nuc*n_gridpoints + lower - 1 + 1]
#                 low_idx = nuc*n_gridpoints + lower - 1 + 1
#             else
#                 low = nuclide_grids[nuc*n_gridpoints + lower + 1]
#                 low_idx = nuc*n_gridpoints + lower + 1
#             end
#         end

#     high = nuclide_grids[low_idx + 1]
#     # Original: high = nuclide_grids[nuc*n_gridpoints + (findfirst(==(low), nuclide_grids) + 1)]
    
#     # calculate the re-useable interpolation factor

#     f = (high.energy - p_energy) / (high.energy - low.energy)
#     if lookups == 2 && j == 4
#         a[37] = f
#     end
#     a[9] = a[9] + 1
    
#     a[16] = high.energy - p_energy
#     a[17] = high.energy - low.energy

#     if lookups == 1
#         a[j + 63] = high.total_xs - f * (high.total_xs - low.total_xs)
#         if j == 5
#             a[95] += a[j + 63] * conc
#             a[96] = conc
#         end
    
#     end
#     if lookups == 2
#         a[j + 73] = high.total_xs - f * (high.total_xs - low.total_xs)
#         if j == 4
#             a[97] += a[j + 73] * conc
#             a[98] = conc
#         end
#     end

#     # Total XS
#     xs_vector[1] = high.total_xs - f * (high.total_xs - low.total_xs)

#     # Elastic XS
#     xs_vector[2] = high.elastic_xs - f * (high.elastic_xs - low.elastic_xs)

#     # Absorption XS
#     xs_vector[3] = high.absorbtion_xs - f * (high.absorbtion_xs - low.absorbtion_xs)

#     # Fission XS
#     xs_vector[4] = high.fission_xs - f * (high.fission_xs - low.fission_xs)

#     # Nu Fission XS
#     xs_vector[5] = high.nu_fission_xs - f * (high.nu_fission_xs - low.nu_fission_xs)


# end




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
    if grid_type == 0 #OriginallY NUCLIDE
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


function add()
    a = 1 + 1
end





end # module

