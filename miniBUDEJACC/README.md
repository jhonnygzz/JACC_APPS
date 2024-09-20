Steps to run the code.

Make sure you are in the miniBUDEJACC directory.

Run these commands: 
```shell
julia --project=. -e 'import Pkg; Pkg.instantiate()'

julia --project=.
include("src/BasicBUDEPreferences.jl")
using .BasicBUDEPreferences

BasicBUDEPreferences.set_backend("cuda")
println(BasicBUDEPreferences.get_backend())

exit()

julia --project=. -e 'include("src/JACCBUDE.jl")'
or
julia --project=. src/JACCBUDE.jl
```

Run profiler ([CUDA](https://cuda.juliagpu.org/stable/development/profiling/)):
```shell
nsys profile julia --project=. -e 'include("src/JACCBUDE.jl")'
nsys stats report1.nsys-rep
```

