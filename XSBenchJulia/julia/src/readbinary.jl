using Base: read, open


struct NuclideGridPoint
    energy::Float64
    total_xs::Float64
    elastic_xs::Float64
    absorbtion_xs::Float64
    fission_xs::Float64
    nu_fission_xs::Float64
end

typedef struct{
	int * num_nucs;                     // Length = length_num_nucs;
	double * concs;                     // Length = length_concs
	int * mats;                         // Length = length_mats
	double * unionized_energy_array;    // Length = length_unionized_energy_array
	int * index_grid;                   // Length = length_index_grid
	NuclideGridPoint * nuclide_grid;    // Length = length_nuclide_grid
    double * p_energy_samples;
    int * mat_samples;

	int length_num_nucs;
	int length_concs;
	int length_mats;
	int length_unionized_energy_array;
	long length_index_grid;
	int length_nuclide_grid;
	int length_p_energy_samples;
	int length_mat_samples;
    int max_num_nucs;

} SimulationData;

function binary_read()
    # println("Current working directory: ", pwd())
    fname = "/home/jvalglz/Fall2024/XSBench/julia/src/XS_data.dat"
    println("Reading all data structures from binary file $fname...")

    open(fname, "r") do fp
        # Read SimulationData Object
        SD = SimulationData(
            read(fp, Int),  # length_num_nucs
            Vector{Int}(undef, read(fp, Int)),  # num_nucs
            read(fp, Int),  # length_concs
            Vector{Float64}(undef, read(fp, Int)),  # concs
            read(fp, Int),  # length_mats
            Vector{Int}(undef, read(fp, Int)),  # mats
            read(fp, Int),  # length_nuclide_grid
            Vector{NuclideGridPoint}(undef, read(fp, Int)),  # nuclide_grid
            read(fp, Int),  # length_index_grid
            Vector{Int}(undef, read(fp, Int)),  # index_grid
            read(fp, Int),  # length_unionized_energy_array
            Vector{Float64}(undef, read(fp, Int))  # unionized_energy_array
        )

        # Read heap arrays into SimulationData Object
        read!(fp, SD.num_nucs)
        read!(fp, SD.concs)
        read!(fp, SD.mats)
        read!(fp, SD.nuclide_grid)
        read!(fp, SD.index_grid)
        read!(fp, SD.unionized_energy_array)
    end

    return SD
end



function print_file_contents(file_path::String)
    try
        # Open the file and read its contents
        open(file_path, "r") do file
            contents = read(file, String)
            println(contents)
        end
    catch ex
        # Handle any errors that may occur
        println("Error: $ex")
    end
end


# Example usage
function main()
    SD = binary_read()
    print_file_contents("/home/jvalglz/Fall2024/XSBench/julia/src/hello.txt")

    # Use the SimulationData object here
    # ...
end

main()
