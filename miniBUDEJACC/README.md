Steps to run the code.

Make sure you are in the miniBUDEJACC directory.

Run these commands: 
```shell
julia --project=. -e 'import Pkg; Pkg.instantiate()'

julia --project=.
using JACC.JACCPreferences
JACCPreferences.set_backend("amdgpu")
exit()

julia --project=.
using JACC
@show JACC.JACCPreferences.backend # should show "amdgpu"
exit()

julia --project=. src/JACCBUDE.jl -p 2 --deck src/data/bm2
```

Run profiler ([CUDA](https://cuda.juliagpu.org/stable/development/profiling/)):
```shell
nsys profile julia --project=. src/JACCBUDE.jl -p 2 --deck src/data/bm2
nsys stats report1.nsys-rep
```

