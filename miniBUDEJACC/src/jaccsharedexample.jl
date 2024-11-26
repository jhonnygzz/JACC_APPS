# import JACC
# import Pkg
# using Test
# using StaticArrays
# include("BUDE.jl")

# @static if endswith(BasicBUDEPreferences.backend, "cuda")
#     # @TODO Julia Pkg.add will add target = :weakdeps in later versions
#     Pkg.add(; name = "CUDA", version = "v5.4.3")

#     #Pkg.add("CUDA")
#     #Pkg.update("CUDA")
#     import CUDA
#     using CUDA # to use cuprintln
#     println("Using CUDA as back end")
#     println("CUDA version: ", CUDA.versioninfo())

#     device = CUDA.device()
#     println("device: $(CUDA.name(device))")

# elseif endswith(BasicBUDEPreferences.backend, "amdgpu")
#     #Pkg.add(; name = "AMDGPU", version = "v0.8.6")
#     Pkg.add("AMDGPU")
#     #Pkg.update("AMDGPU")
#     import AMDGPU
#     println("Using AMDGPU as back end")
    
#     device = AMDGPU.device()
#     println("device: $device")

#   elseif endswith(BasicBUDEPreferences.backend, "threads")
#     println("Using threads as back end")
# end


# # Define init_image function
# function init_image(type, num_bands, num_voxel, size_voxel)
#     return rand(type, num_bands, num_voxel, size_voxel)
# end

# # Define init_filter function
# function init_filter(type, size_voxel)
#     return rand(type, size_voxel)
# end

# # Define spectral_shared function
# function spectral_shared(i, j, image, filter, num_bands)
#     # Shared memory initialization
#     filter_shared = JACC.shared(filter)
#     for b in 1:num_bands
#         @inbounds image[b, i, j] *= filter_shared[j]
#     end
# end

# # Main code
# num_bands = 60
# num_voxel = 10_240
# size_voxel = 64*64

# image = init_image(Float32, num_bands, num_voxel, size_voxel)
# filter = init_filter(Float32, size_voxel)

# jimage = JACC.Array(image)
# jfilter = JACC.Array(filter)

# JACC.parallel_for((num_voxel,size_voxel), spectral_shared, jimage, jfilter, num_bands)