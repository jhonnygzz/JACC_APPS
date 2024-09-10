Steps to run the code.

Make sure you are in the miniBUDEJACC directory.

Run these commands: 
```shell
julia --project=. -e 'import Pkg; Pkg.instantiate()'

julia --project=.
include("src/BasicBUDEPreferences.jl")
using .BasicBUDEPreferences

BasicBUDEPreferences.set_backend("jacc-cuda")
println(BasicBUDEPreferences.get_backend())

exit()

julia --project=. -e 'include("src/JACCBUDE.jl")'
```
