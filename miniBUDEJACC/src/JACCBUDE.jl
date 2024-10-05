# IMPORTANT NOTES

#TODO: Update github repository on how to run the code: DONE
#Compare apples to apples
#Start with the threads version always for code translation to JACC
#TODO: Check its using elapsed correctly: DONE
#TODO: Print the time: DONE
#TODO: Make sure the code is using the CUDA device: DONE
#TODO: First compare the original versions CUDA vs Threads vs AMDGPU: Done
# We cannot have different if statements (only one function for the kernel)
# Play with mi250x in the future
# Everything related to calculations uses Float32 (Single Precision). For output related, it uses Float64 (Double precision), examples include Kernel time, Average time, etc.
#TODO: Should the time calculation be exactly before the parallel_for?

# There is a process in execution when running this code in CUDA (A100):
# miniBUDE CUDA: 0   N/A  N/A   1374566      C   julia                                         446MiB
# miniBUDEJACC CUDA:  0   N/A  N/A   1378053      C   julia                                        1134MiB
# Ran a profiler (NVIDIA Nsight Systems)
# Had to do CUDA.synchronize() to get correct times when using CUDA. This does not happen with AMDGPU
# Difference between the original CUDA code and this code may be the work distribution affecting the order of operations
# Can't do print statements of variables inside the kernel function since it is running on the GPU
# Already compared the output of both JACC-CUDA and CUDA. The output is mostly the same, expect some energies are 0.01 different.
# TODO: Print the values of the variables. (protein, forcefield, etc)
# after updating jacc to jacc#main, now the it synchronizes correctly without me having to do cuda.synchronize()



import JACC
import Pkg
using Test
using StaticArrays
#using JACC
include("BUDE.jl")


include("BasicBUDEPreferences.jl")
@static if endswith(BasicBUDEPreferences.backend, "cuda")
    # @TODO Julia Pkg.add will add target = :weakdeps in later versions
    Pkg.add(; name = "CUDA", version = "v5.4.3")

    #Pkg.add("CUDA")
    #Pkg.update("CUDA")
    import CUDA
    using CUDA # to use cuprintln
    println("Using CUDA as back end")
    println("CUDA version: ", CUDA.versioninfo())

    device = CUDA.device()
    println("device: $(CUDA.name(device))")

elseif endswith(BasicBUDEPreferences.backend, "amdgpu")
    #Pkg.add(; name = "AMDGPU", version = "v0.8.6")
    Pkg.add("AMDGPU")
    #Pkg.update("AMDGPU")
    import AMDGPU
    println("Using AMDGPU as back end")
    
    device = AMDGPU.device()
    println("device: $device")

  elseif endswith(BasicBUDEPreferences.backend, "threads")
    println("Using threads as back end")
end

# const Device = (undef, "CPU", "Threaded")

# function devices()
#   return [Device]
# end



const Float32Max = typemax(Float32)

function run(params::Params, deck::Deck) #_::DeviceWithRepr)
  #println("Using max $(Threads.nthreads()) threads")
  # if params.ppwi != DefaultPPWI
  #   @warn "Threaded implementation only uses wgsize, the PPWI argument is ignored"
  # end

  #nposes::Int = size(deck.poses)[2]
  #poses = size(deck.poses)[2]
  protein = JACC.Array{Atom}(deck.protein) #The number of atoms in the protein molecule. The protein is typically the larger molecule that the ligand is being docked to.
  ligand = JACC.Array{Atom}(deck.ligand) #The number of atoms in the ligand molecule. The ligand is typically a small molecule that we're trying to dock to the prot
  forcefield = JACC.Array{FFParams}(deck.forcefield) # The number of force field parameters used in the simulation. Force fields define the interactions between atoms and are crucial for calculating the energy of each pose.
  poses = JACC.Array{Float32,2}(deck.poses) #This indicates the number of different orientations and positions of the ligand molecule being tested
  #etotals = JACC.Array{Float32}(undef, poses)

  #forcefield = JACC.shared(forcefield)
  #println("worked")

  etotals = JACC.Array{Float32}(undef, size(deck.poses)[2]) #stores the total energy for each pose.

  #CUDA.@profile sin.(etotals)

  #println("Type of protein: ", typeof(protein))



  # warmup
  fasten_main(
    Val(convert(Int, params.wgsize)),
    Val(convert(Int, params.ppwi)),
    protein,
    ligand,
    forcefield,
    poses,
    etotals,
  )


  # ORIGINAL ELAPSED FROM THREADED.JL

  # elapsed = @elapsed for _ = 1:params.iterations
  #   fasten_main(
  #     Val(convert(Int, params.wgsize)),
  #     protein,
  #     ligand,
  #     forcefield,
  #     poses,
  #     etotals,
  #   )
  # end


  #MODIFIED ELAPSED FROM THREADED.JL
  total_elapsed = 0.0
  for i = 1:params.iterations
    start_time = time()
    begin
      fasten_main(
        Val(convert(Int, params.wgsize)),
        Val(convert(Int, params.ppwi)),
        protein,
        ligand,
        forcefield,
        poses,
        etotals,
      )
    end
    #CUDA.synchronize()
    #AMDGPU.synchronize()
    end_time = time()
    iteration_elapsed = end_time - start_time
    total_elapsed += iteration_elapsed
    println("Iteration $(i): $(iteration_elapsed) seconds")
  end

  average_time = total_elapsed / params.iterations
  println("Total elapsed time: $(total_elapsed) seconds")
  println("Average time per iteration: $(average_time) seconds")
  

  # Transfer etotals back to CPU (for cuda and amdgpu only)
  etotals_cpu = Array(etotals)

  
  #(etotals, total_elapsed, params.wgsize)
  #(etotals_cpu, elapsed, params.wgsize)
  (etotals_cpu, total_elapsed, params.wgsize)
end


@fastmath function fasten_main(
  ::Val{WGSIZE},
  ::Val{PPWI},
  protein::JACC.Array{Atom},
  ligand::JACC.Array{Atom},
  forcefield::JACC.Array{FFParams},
  poses::JACC.Array{Float32,2},
  etotals::JACC.Array{Float32},
) where {WGSIZE, PPWI}

  nposes::Int = size(poses)[2]
  nligand::Int = length(ligand)
  nprotein::Int = length(protein)
  total_work_items::Int = ceil(Int, nposes / PPWI)
  #print("blocks_size:", blocks_size)
  

  function kernel(index)
    # Get index of first TD (Thread Data)
    ix = (index - 1) * PPWI + 1 #alculates the starting pose index for each work item
    ix = ix <= nposes ? ix : nposes - PPWI #make sure we dont go out of bounds (not go past)

    etot = MArray{Tuple{PPWI},Float32}(undef)
    transform = MArray{Tuple{PPWI,3},Vec4f32}(undef)

    @simd for i = 1:PPWI
      pose_index = ix + i - 1
      @inbounds sx::Float32 = sin(poses[1, pose_index])
      @inbounds cx::Float32 = cos(poses[1, pose_index])
      @inbounds sy::Float32 = sin(poses[2, pose_index])
      @inbounds cy::Float32 = cos(poses[2, pose_index])
      @inbounds sz::Float32 = sin(poses[3, pose_index])
      @inbounds cz::Float32 = cos(poses[3, pose_index])
      @inbounds transform[i, 1] = #
        Vec4f32(cy * cz, sx * sy * cz - cx * sz, cx * sy * cz + sx * sz, poses[4, pose_index])
      @inbounds transform[i, 2] = #
        Vec4f32(cy * sz, sx * sy * sz + cx * cz, cx * sy * sz - sx * cz, poses[5, pose_index])
      @inbounds transform[i, 3] = #
        Vec4f32(-sy, sx * cy, cx * cy, poses[6, pose_index])
      @inbounds etot[i] = 0
    end

    #for il = 1:nligand
    @inbounds for l_atom::Atom in ligand
      #l_atom = ligand[il]
      l_params = forcefield[l_atom.type+1]
      lhphb_ltz = l_params.hphb < Zero
      lhphb_gtz = l_params.hphb > Zero

      lpos = MArray{Tuple{PPWI},Vec3f32}(undef)

    @simd for i = 1:PPWI
      # Transform ligand atom
      @inbounds lpos[i] = Vec3f32(
        transform[i, 1].w +
        l_atom.x * transform[i, 1].x +
        l_atom.y * transform[i, 1].y +
        l_atom.z * transform[i, 1].z,
        transform[i, 2].w +
        l_atom.x * transform[i, 2].x +
        l_atom.y * transform[i, 2].y +
        l_atom.z * transform[i, 2].z,
        transform[i, 3].w +
        l_atom.x * transform[i, 3].x +
        l_atom.y * transform[i, 3].y +
        l_atom.z * transform[i, 3].z,
      )
    end

      #for ip = 1:nprotein
      @inbounds for p_atom::Atom in protein
        #p_atom = protein[ip]
        #p_params = forcefield[p_atom.type+1]
        @inbounds p_params::FFParams = forcefield[p_atom.type+1]


        radij::Float32 = p_params.radius + l_params.radius
        r_radij::Float32 = One / radij
  
        elcdst::Float32 = (p_params.hbtype == HbtypeF && l_params.hbtype == HbtypeF) ? Four : Two
        elcdst1::Float32 = (p_params.hbtype == HbtypeF && l_params.hbtype == HbtypeF) ? Quarter : Half
        type_E::Bool = ((p_params.hbtype == HbtypeE || l_params.hbtype == HbtypeE))
  
        phphb_ltz::Bool = p_params.hphb < Zero
        phphb_gtz::Bool = p_params.hphb > Zero
        phphb_nz::Bool = p_params.hphb != Zero
        p_hphb::Float32 = p_params.hphb * (phphb_ltz && lhphb_gtz ? -One : One)
        l_hphb::Float32 = l_params.hphb * (phphb_gtz && lhphb_ltz ? -One : One)
        distdslv::Float32 =
          (phphb_ltz ? (lhphb_ltz ? Npnpdist : Nppdist) : (lhphb_ltz ? Nppdist : -Float32Max))
        r_distdslv::Float32 = One / distdslv
  
        chrg_init::Float32 = l_params.elsc * p_params.elsc
        dslv_init::Float32 = p_hphb + l_hphb

      @simd for i = 1:PPWI
        @inbounds x::Float32 = lpos[i].x - p_atom.x
        @inbounds y::Float32 = lpos[i].y - p_atom.y
        @inbounds z::Float32 = lpos[i].z - p_atom.z

        distij::Float32 = sqrt(x ^ 2 + y ^ 2 + z ^ 2)

        # Calculate the sum of the sphere radii
        distbb::Float32 = distij - radij
        zone1::Bool = (distbb < Zero)

        # Calculate steric energy
        @inbounds etot[i] += (One - (distij * r_radij)) * (zone1 ? Two * Hardness : Zero)

        # Calculate formal and dipole charge interactions
        chrg_e::Float32 =
          chrg_init * ((zone1 ? One : (One - distbb * elcdst1)) * (distbb < elcdst ? One : Zero))
        neg_chrg_e::Float32 = -abs(chrg_e)
        chrg_e = type_E ? neg_chrg_e : chrg_e
        @inbounds etot[i] += chrg_e * Cnstnt

        # Calculate the two cases for Nonpolar-Polar repulsive interactions
        coeff::Float32 = (One - (distbb * r_distdslv))
        dslv_e::Float32 = dslv_init * ((distbb < distdslv && phphb_nz) ? One : Zero)
        dslv_e *= (zone1 ? One : coeff)
        @inbounds etot[i] += dslv_e
      end
    end
  end

    @simd for i = 1:PPWI
      pose_index = ix + i - 1
      if pose_index <= nposes
        @inbounds etotals[pose_index] = etot[i] * Half
      end
    end
  end #end for kernel function

  JACC.parallel_for(total_work_items, kernel)
end







#   function kernel(index, protein, ligand, forcefield, poses, etotals)
#     #ix = (index - 1) + 1
#     ix = (index - 1) * PPWI + 1

#     ix = ix <= nposes ? ix : nposes - PPWI #out-of-bounds access (ensures that no thread tries to access poses beyond the total number of poses (nposes))
#     #ix = min(ix, nposes - PPWI + 1)

#     etot = MArray{Tuple{PPWI}, Float32}(undef)
#     transform = MArray{Tuple{PPWI, 3, 4}, Float32}(undef)

#     for i = 1:PPWI
#       pose_index = ix + i - 1
#       sx = sin(poses[1, pose_index])
#       cx = cos(poses[1, pose_index])
#       sy = sin(poses[2, pose_index])
#       cy = cos(poses[2, pose_index])
#       sz = sin(poses[3, pose_index])
#       cz = cos(poses[3, pose_index])

#       transform[i, 1, 1] = cy * cz
#       transform[i, 1, 2] = sx * sy * cz - cx * sz
#       transform[i, 1, 3] = cx * sy * cz + sx * sz
#       transform[i, 1, 4] = poses[4, pose_index]
#       transform[i, 2, 1] = cy * sz
#       transform[i, 2, 2] = sx * sy * sz + cx * cz
#       transform[i, 2, 3] = cx * sy * sz - sx * cz
#       transform[i, 2, 4] = poses[5, pose_index]
#       transform[i, 3, 1] = -sy
#       transform[i, 3, 2] = sx * cy
#       transform[i, 3, 3] = cx * cy
#       transform[i, 3, 4] = poses[6, pose_index]
#       etot[i] = Zero
#     end

#     for l_atom in ligand
#       l_params = forcefield[l_atom.type+1]
#       lhphb_ltz = l_params.hphb < Zero
#       lhphb_gtz = l_params.hphb > Zero

#       lpos = MArray{Tuple{PPWI, 3}, Float32}(undef)

#       for i = 1:PPWI
#         lpos[i, 1] = (
#           transform[i, 1, 4] +
#           l_atom.x * transform[i, 1, 1] +
#           l_atom.y * transform[i, 1, 2] +
#           l_atom.z * transform[i, 1, 3]
#         )
#         lpos[i, 2] = (
#           transform[i, 2, 4] +
#           l_atom.x * transform[i, 2, 1] +
#           l_atom.y * transform[i, 2, 2] +
#           l_atom.z * transform[i, 2, 3]
#         )
#         lpos[i, 3] = (
#           transform[i, 3, 4] +
#           l_atom.x * transform[i, 3, 1] +
#           l_atom.y * transform[i, 3, 2] +
#           l_atom.z * transform[i, 3, 3]
#         )
#       end

#       for p_atom in protein
#         p_params = forcefield[p_atom.type+1]

#         radij = p_params.radius + l_params.radius
#         r_radij = One / radij

#         elcdst = (p_params.hbtype == HbtypeF && l_params.hbtype == HbtypeF) ? Four : Two
#         elcdst1 = (p_params.hbtype == HbtypeF && l_params.hbtype == HbtypeF) ? Quarter : Half
#         type_E = ((p_params.hbtype == HbtypeE || l_params.hbtype == HbtypeE))

#         phphb_ltz = p_params.hphb < Zero
#         phphb_gtz = p_params.hphb > Zero
#         phphb_nz = p_params.hphb != Zero
#         p_hphb = p_params.hphb * (phphb_ltz && lhphb_gtz ? -One : One)
#         l_hphb = l_params.hphb * (phphb_gtz && lhphb_ltz ? -One : One)
#         distdslv = (phphb_ltz ? (lhphb_ltz ? Npnpdist : Nppdist) : (lhphb_ltz ? Nppdist : -Float32Max))
#         r_distdslv = One / distdslv

#         chrg_init = l_params.elsc * p_params.elsc
#         dslv_init = p_hphb + l_hphb

#         for i = 1:PPWI
#           x = lpos[i, 1] - p_atom.x
#           y = lpos[i, 2] - p_atom.y
#           z = lpos[i, 3] - p_atom.z

#           distij = sqrt(x * x + y * y + z * z)

#           distbb = distij - radij
#           zone1 = (distbb < Zero)

#           etot[i] += (One - (distij * r_radij)) * (zone1 ? Two * Hardness : Zero)

#           chrg_e = chrg_init * ((zone1 ? One : (One - distbb * elcdst1)) * (distbb < elcdst ? One : Zero))
#           neg_chrg_e = -abs(chrg_e)
#           chrg_e = type_E ? neg_chrg_e : chrg_e
#           etot[i] += chrg_e * Cnstnt

#           coeff = (One - (distbb * r_distdslv))
#           dslv_e = dslv_init * ((distbb < distdslv && phphb_nz) ? One : Zero)
#           dslv_e *= (zone1 ? One : coeff)
#           etot[i] += dslv_e
#         end
#       end
#     end

#     for i = 1:PPWI
#       pose_index = ix + i - 1
#       if pose_index <= nposes
#         etotals[pose_index] = etot[i] * Half
#       end
#     end
#   end

#   JACC.parallel_for(blocks_size, kernel, protein, ligand, forcefield, poses, etotals)
# end

























# @fastmath function fasten_main(
#   ::Val{WGSIZE},
#   protein::JACC.Array{Atom},
#   ligand::JACC.Array{Atom},
#   forcefield::JACC.Array{FFParams},
#   poses::JACC.Array{Float32,2},
#   etotals::JACC.Array{Float32},
# ) where {WGSIZE} #wgsize can be any type, but it must be specified when the function is called 
# # can be called like this fasten_main(Val(64), protein, ligand, forcefield, poses, etotals) which means work group size of 64.

#   #local variables
#   nposes::Int = size(poses)[2]
#   numGroups::Int = nposes ÷ WGSIZE
#   nligand::Int = length(ligand)
#   nprotein::Int = length(protein)

#   #println("Printing value of forcefield: ", forcefield)

#   #println("Type of forcefield: ", typeof(forcefield))


#   # if isdefined(JACC, :shared)
#   #   println("JACC.shared exists")
#   # else
#   #   println("JACC.shared does not exist")
#   # end
#   #counter = JACC.Array{Int32}(undef, 1)
#   #JACC.fill!(counter, 0)  # Initialize the counter to 0


#   #Threads.@threads for group = 1:numGroups (use a function instead)

#   function kernel(group, protein, ligand, forcefield, poses, etotals) #counter)

#     #utilizing shared memory
#     #JACC.shared(ligand) #this also works, which does not make sense.
#     #forcefield_shared = JACC.shared(forcefield)
#     #JACC.@atomic counter[1] += 1
  
  

#     etot = MArray{Tuple{WGSIZE}, Float32}(undef)
#     transform = MArray{Tuple{WGSIZE, 3, 4},Float32}(undef)

#     @simd for i = 1:WGSIZE
#       ix = (group - 1) * (WGSIZE) + i
#       @inbounds sx::Float32 = sin(poses[1, ix])
#       @inbounds cx::Float32 = cos(poses[1, ix])
#       @inbounds sy::Float32 = sin(poses[2, ix])
#       @inbounds cy::Float32 = cos(poses[2, ix])
#       @inbounds sz::Float32 = sin(poses[3, ix])
#       @inbounds cz::Float32 = cos(poses[3, ix])

#       # println("Iteration $i:")
#       # println("  ix: $ix")
#       # println("  sx: $sx, cx: $cx")
#       # println("  sy: $sy, cy: $cy")
#       # println("  sz: $sz, cz: $cz")
#       #@cuprintln("Group: ", group)
#       @inbounds transform[i, 1, 1] = cy * cz
#       @inbounds transform[i, 1, 2] = sx * sy * cz - cx * sz
#       @inbounds transform[i, 1, 3] = cx * sy * cz + sx * sz
#       @inbounds transform[i, 1, 4] = poses[4, ix]
#       @inbounds transform[i, 2, 1] = cy * sz
#       @inbounds transform[i, 2, 2] = sx * sy * sz + cx * cz
#       @inbounds transform[i, 2, 3] = cx * sy * sz - sx * cz
#       @inbounds transform[i, 2, 4] = poses[5, ix]
#       @inbounds transform[i, 3, 1] = -sy
#       @inbounds transform[i, 3, 2] = sx * cy
#       @inbounds transform[i, 3, 3] = cx * cy
#       @inbounds transform[i, 3, 4] = poses[6, ix]
#       @inbounds etot[i] = Zero
#     end

#     for il = 1:nligand
#       @inbounds l_atom::Atom = ligand[il]

#       @inbounds l_params::FFParams = forcefield[l_atom.type+1]
#       lhphb_ltz::Bool = l_params.hphb < Zero
#       lhphb_gtz::Bool = l_params.hphb > Zero

#       lpos = MArray{Tuple{WGSIZE, 3}, Float32}(undef)

#       @simd for i = 1:WGSIZE
#         @inbounds lpos[i, 1] = (
#           transform[i, 1, 4] +
#           l_atom.x * transform[i, 1, 1] +
#           l_atom.y * transform[i, 1, 2] +
#           l_atom.z * transform[i, 1, 3]
#         )
#         @inbounds lpos[i, 2] = (
#           transform[i, 2, 4] +
#           l_atom.x * transform[i, 2, 1] +
#           l_atom.y * transform[i, 2, 2] +
#           l_atom.z * transform[i, 2, 3]
#         )
#         @inbounds lpos[i, 3] = (
#           transform[i, 3, 4] +
#           l_atom.x * transform[i, 3, 1] +
#           l_atom.y * transform[i, 3, 2] +
#           l_atom.z * transform[i, 3, 3]
#         )
#       end

#       for ip = 1:nprotein
#         @inbounds p_atom = protein[ip]
#         @inbounds p_params::FFParams = forcefield[p_atom.type+1]

#         radij::Float32 = p_params.radius + l_params.radius
#         r_radij::Float32 = One / radij

#         elcdst::Float32 = (p_params.hbtype == HbtypeF && l_params.hbtype == HbtypeF) ? Four : Two
#         elcdst1::Float32 =
#           (p_params.hbtype == HbtypeF && l_params.hbtype == HbtypeF) ? Quarter : Half
#         type_E::Bool = ((p_params.hbtype == HbtypeE || l_params.hbtype == HbtypeE))

#         phphb_ltz::Bool = p_params.hphb < Zero
#         phphb_gtz::Bool = p_params.hphb > Zero
#         phphb_nz::Bool = p_params.hphb != Zero
#         p_hphb::Float32 = p_params.hphb * (phphb_ltz && lhphb_gtz ? -One : One)
#         l_hphb::Float32 = l_params.hphb * (phphb_gtz && lhphb_ltz ? -One : One)
#         distdslv::Float32 =
#           (phphb_ltz ? (lhphb_ltz ? Npnpdist : Nppdist) : (lhphb_ltz ? Nppdist : -Float32Max))
#         r_distdslv::Float32 = One / distdslv

#         chrg_init::Float32 = l_params.elsc * p_params.elsc
#         dslv_init::Float32 = p_hphb + l_hphb

#         @simd for i = 1:WGSIZE
#           @inbounds x::Float32 = lpos[i, 1] - p_atom.x
#           @inbounds y::Float32 = lpos[i, 2] - p_atom.y
#           @inbounds z::Float32 = lpos[i, 3] - p_atom.z

#           distij::Float32 = sqrt(x ^ 2 + y ^ 2+ z ^ 2)

#           # Calculate the sum of the sphere radii
#           distbb::Float32 = distij - radij
#           zone1::Bool = (distbb < Zero)

#           # Calculate steric energy
#           @inbounds etot[i] += (One - (distij * r_radij)) * (zone1 ? Two * Hardness : Zero)

#           # Calculate formal and dipole charge interactions
#           chrg_e::Float32 =
#             chrg_init * ((zone1 ? One : (One - distbb * elcdst1)) * (distbb < elcdst ? One : Zero))
#           neg_chrg_e::Float32 = -abs(chrg_e)
#           chrg_e = type_E ? neg_chrg_e : chrg_e
#           @inbounds etot[i] += chrg_e * Cnstnt

#           # Calculate the two cases for Nonpolar-Polar repulsive interactions
#           coeff::Float32 = (One - (distbb * r_distdslv))
#           dslv_e::Float32 = dslv_init * ((distbb < distdslv && phphb_nz) ? One : Zero)
#           dslv_e *= (zone1 ? One : coeff)
#           @inbounds etot[i] += dslv_e
#         end
#       end
#     end
#     @simd for i = 1:WGSIZE
#       ix = (group - 1) * WGSIZE + i
#       @inbounds etotals[ix] = etot[i] * Half
#     end
#   end

#   # total_elapsed = 0.0
#   # start_time = time()
#   #println("Kernel will be called $(numGroups) times")
#   JACC.parallel_for(numGroups, kernel, protein, ligand, forcefield, poses, etotals) #counter) #numgroups (first variable) determines how many times the kernel function will be called in parallel
#   #host_counter = Array(counter)
#   #println("Kernel was actually called $(host_counter[1]) times")
#   #println("Kernel was actually called $(counter) times")
#   # CUDA.synchronize()
#   # end_time = time()
#   # iteration_elapsed = end_time - start_time
#   # total_elapsed += iteration_elapsed
#   # println("Iteration parallel_for: $(iteration_elapsed) seconds")

# end




main()
