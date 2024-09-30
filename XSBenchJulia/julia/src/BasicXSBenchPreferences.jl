module BasicXSBenchPreferences

using Preferences

const PACKAGE_UUID = Base.UUID("a49ee08d-2cae-4ecc-ace2-0baa55f7875c")  # Replace with your generated UUID

# Function to set the backend preference
function set_backend(new_backend::String)
    new_backend_lc = lowercase(new_backend)
    if !(new_backend_lc in ("jacc-threads", "jacc-cuda", "jacc-amdgpu"))
        throw(ArgumentError("Invalid backend: \"$(new_backend)\" accepted values: jacc-threads, jacc-cuda, jacc-amdgpu"))
    end

    # Set the preference
    set_preferences!(PACKAGE_UUID, "backend" => new_backend_lc, force=true)
    @info("New backend set; restart your Julia session for this change to take effect!")
end

# Function to get the backend preference
function get_backend()
    Preferences.load_preference(PACKAGE_UUID, "backend", "amdgpu")
end

# Load the backend preference
const backend = get_backend()

# Print the current backend preference
println("backend: ", backend)

end # module BasicXSBenchPreferences