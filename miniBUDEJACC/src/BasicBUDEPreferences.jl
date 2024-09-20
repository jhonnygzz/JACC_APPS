module BasicBUDEPreferences

using Preferences

const PACKAGE_UUID = Base.UUID("a49ee08d-2cae-4ecc-ace2-0baa55f7875c")

# taken from https://github.com/JuliaPackaging/Preferences.jl
function set_backend(new_backend::String)
    new_backend_lc = lowercase(new_backend)
    if !(new_backend_lc in ("threads", "cuda", "amdgpu"))
        throw(ArgumentError("Invalid backend: \"$(new_backend)\" accepted values: threads, cuda, amdgpu"))
    end
    # Set it in our runtime values, as well as saving it to disk
    set_preferences!(PACKAGE_UUID, "backend" => new_backend_lc, force=true)
    #@set_preferences!("backend"=>new_backend_lc)
    @info("New backend set; restart your Julia session for this change to take effect!")
end

function get_backend()
   Preferences.load_preference(PACKAGE_UUID, "backend", "threads")
end

const backend = get_backend()
#const backend = @load_preference("backend", "threads")

#export set_backend, get_backend, backend

end # module BasicBUDEPreferences