module BasicBUDEPreferences

using Preferences

const PACKAGE_UUID = Base.UUID("a49ee08d-2cae-4ecc-ace2-0baa55f7875c")

function set_backend(new_backend::String)
    new_backend_lc = lowercase(new_backend)
    if !(new_backend_lc in ("jacc-threads", "jacc-cuda", "jacc-amdgpu"))
        throw(ArgumentError("Invalid backend: \"$(new_backend)\" accepted values: jacc-threads, jacc-cuda, jacc-amdgpu"))
    end
    # Set it in our runtime values, as well as saving it to disk
    set_preferences!(PACKAGE_UUID, "backend" => new_backend_lc, force=true)
    @info("New backend set; restart your Julia session for this change to take effect!")
end

function get_backend()
    Preferences.load_preference(PACKAGE_UUID, "backend", "jacc-threads")
end

const backend = get_backend()

#export set_backend, get_backend, backend

end # module BasicBUDEPreferences