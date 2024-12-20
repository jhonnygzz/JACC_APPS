include("Stream.jl")
import JACC
import Pkg
using JACC
# const JACCData = StreamData{T,JACC.Array{T}} where {T}

const TBSize = 1024::Int
const DotBlocks = 256::Int

# set backend
include("BasicBabelPreferences.jl")
backend = "threads"

@static if endswith(BasicBabelPreferences.backend, "cuda")
    # @TODO Julia Pkg.add will add target = :weakdeps in later versions
    Pkg.add(; name = "CUDA", version = "v5.4.3")
    
    # Pkg.add("CUDA")
    # Pkg.update("CUDA")
    import CUDA
    # using CUDA
    CUDA.set_runtime_version!(v"12.1.0")  # Replace with a supported CUDA version

    # pkg = Base.PkgId(Base.UUID("76a88914-d11a-5bdc-97e0-2f5a05c973a2"), "CUDA_Runtime_jll")
    # Base.compilecache(pkg)

    # CUDA.set_runtime_version!(v"11.2")  # Replace with your installed CUDA version
    # println("Using CUDA as back end")
    # println("CUDA version: ", CUDA.versioninfo())

    device = CUDA.device()
    # println("device: $(CUDA.name(device))")
    backend = "cuda"

elseif endswith(BasicBabelPreferences.backend, "amdgpu")
    # Pkg.add(; name = "AMDGPU", version = "v0.8.11")
    # Pkg.add(; name = "AMDGPU", version = "v1.0.3") 
    import AMDGPU
    Pkg.add("AMDGPU")
    Pkg.update("AMDGPU")
    
    # using AMDGPU
    println("Using AMDGPU as back end")
    
    device = AMDGPU.device()
    println("device: $device")

  elseif endswith(BasicBabelPreferences.backend, "threads")
    using Base.Threads
    println("Using threads as back end")
end

const JACCData = StreamData{T,JACC.Array{T}} where {T}

function devices()::Vector{DeviceWithRepr}
    if backend == "cuda"
        return !CUDA.functional(false) ? String[] :
        map(d -> (d, "$(CUDA.name(d)) ($(repr(d)))", "CUDA.jl"), CUDA.devices())
    elseif backend == "amdgpu"
        try
            # Check if ROCm GPUs are available
            if AMDGPU.has_rocm_gpu()
              println("ROCm GPU detected.")
              
              # Retrieve the list of ROCm devices
              devices = AMDGPU.devices()
              # println("Detected devices: $devices")
              
              # Sort the devices by their string representation
              sorted_devices = sort(devices, by = repr)
              # println("Sorted devices: $sorted_devices")
              
              # Map the sorted devices to the desired format
              device_list = map(x -> (x, repr(x), "AMDGPU.jl"), sorted_devices)
              # println("Device list: $device_list")
              
              return device_list
            else
              println("No ROCm GPU detected.")
              return DeviceWithRepr[]
            end
          catch e
            # Handle any errors that occur during the process
            println("Error in devices(): $e")
            return DeviceWithRepr[]
          end
    else
        return [(undef, "$(Sys.cpu_info()[1].model) ($(Threads.nthreads())T)", "Threaded")]
    end
end
  



function make_stream(
    arraysize::Int,
    scalar::T,
    device::DeviceWithRepr,
    silent::Bool,
)::Tuple{JACCData{T},Nothing} where {T}

    if arraysize % TBSize != 0
        error("arraysize ($(arraysize)) must be divisible by $(TBSize)!")
    end

    if backend == "cuda"
        CUDA.device!(device[1])
        selected = CUDA.device()
        # show_reason is set to true here so it dumps CUDA info 
        # for us regardless of whether it's functional
        if !CUDA.functional(true)
            error("Non-functional CUDA configuration")
        end
        if !silent
            println("Using CUDA device: $(CUDA.name(selected)) ($(repr(selected)))")
            println("Kernel parameters: <<<$(arraysize ÷ TBSize),$(TBSize)>>>")
        end
        
    elseif backend == "amdgpu"
        # XXX AMDGPU doesn't expose an API for setting the default like CUDA.device!()
        # but AMDGPU.get_default_agent returns DEFAULT_AGENT so we can do it by hand
        # AMDGPU.DEFAULT_AGENT[] = device[1]
        # selected = AMDGPU.get_default_agent()
        if !silent
            # println("Using GPU HSA device: $(AMDGPU.get_name(selected)) ($(repr(selected)))")
            # println("Kernel parameters   : <<<$(arraysize),$(TBSize)>>>")
        end
    else
        if !silent
            println("Using CPU device")
            println("Kernel parameters: <<<$(arraysize ÷ TBSize),$(TBSize)>>>")
        end
    end

    # a = JACC.Array{T}(undef, arraysize)
    # b = JACC.Array{T}(undef, arraysize)
    # c = JACC.Array{T}(undef, arraysize)
        return (
        JACCData{T}(
            JACC.Array{T}(undef, arraysize),
            JACC.Array{T}(undef, arraysize),
            JACC.Array{T}(undef, arraysize),
            scalar,
            arraysize,
        ),
        nothing,
    )
    # return data, nothing

end

function init_arrays!(data::JACCData{T}, _, init::Tuple{T,T,T}) where {T}
    JACC.fill!(data.a, init[1])
    JACC.fill!(data.b, init[2])
    JACC.fill!(data.c, init[3])
end

function copy!(data::JACCData{T}, _) where {T}
    # println("Type of data.a in copy! (start): ", typeof(data.a))
    function kernel(i, a::AbstractArray{T}, c::AbstractArray{T})
        @inbounds c[i] = a[i]
        return
    end
    JACC.parallel_for(data.size, kernel, data.a, data.c)
end

function mul!(data::JACCData{T}, _) where {T}
    function kernel(i, b::AbstractArray{T}, c::AbstractArray{T}, scalar::T)
        @inbounds b[i] = scalar * c[i]
        return
    end
    # JACC.parallel_for(data.size, kernel, JACC.Array(data.b), JACC.Array(data.c), data.scalar)
    JACC.parallel_for(data.size, kernel, data.b, data.c, data.scalar)
end

function add!(data::JACCData{T}, _) where {T}
    function kernel(i, a::AbstractArray{T}, b::AbstractArray{T}, c::AbstractArray{T})
        @inbounds c[i] = a[i] + b[i]
        return
    end
    # JACC.parallel_for(data.size, kernel, JACC.Array(data.a), JACC.Array(data.b), JACC.Array(data.c))
    JACC.parallel_for(data.size, kernel, data.a, data.b, data.c)
end

function triad!(data::JACCData{T}, _) where {T}
    function kernel(i, a::AbstractArray{T}, b::AbstractArray{T}, c::AbstractArray{T}, scalar::T)
        @inbounds a[i] = b[i] + (scalar * c[i])
        return
    end
    JACC.parallel_for(data.size, kernel, data.a, data.b,data.c, data.scalar)
end

function dot(data::JACCData{T}, _) where {T}
    function kernel(i, a::AbstractArray{T}, b::AbstractArray{T}, summation::AbstractArray{T})
        # @inbounds c[i] = a[i] * b[i]
        @inbounds summation[i] = a[i] * b[i]
        return
    end

    function reduce(i, summation::AbstractArray{T})
        summation[i]
    end
    summation = JACC.Array{T}(undef, data.size)
    JACC.parallel_for(data.size, kernel, data.a, data.b, summation)
    # summation = JACC.parallel_reduce(data.size, reduce, summation)
    # sum = Base.Array(summation)

    # summation = sum(summation)
    
    return (Base.Array(JACC.parallel_reduce(data.size, reduce, summation)))[1]
end

function read_data(data::JACCData{T}, _)::VectorData{T} where {T}
    return VectorData{T}(data.a, data.b, data.c, data.scalar, data.size)
end



# function free!(data::JACCData{T}) where {T}
#     JACC.free!(data.a)
#     JACC.free!(data.b)
#     JACC.free!(data.c)
# end




main()