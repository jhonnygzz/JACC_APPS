Examples to run the code (JACC, JULIA and C/C++)

JACC:
Make sure you are in the miniBUDEJACC directory.

Run these commands: 
```shell
module load julia
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

C/C++:
Make sure you are in the miniBUDECJulia directory.
```shell
module load nvhpc
nvc++ --version
nvcc --version
cmake -Bbuild_cudaH100 -H. -DMODEL=cuda \
-DCMAKE_CXX_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/24.9/compilers/bin/nvc++ \
-DCMAKE_CUDA_COMPILER=/opt/nvidia/hpc_sdk/Linux_x86_64/24.9/compilers/bin/nvcc \
-DCXX_EXTRA_FLAGS="-march=native -O3 -ffast-math" \
-DCUDA_ARCH=sm_90
cmake --build build_cudaH100
cd build_cudaH100
BEST RUN: ./cuda-bude --deck ../data/bm2 --ppwi 4 --wgsize 128 --device 1
```



JULIA:
Make sure you are in the miniBUDECJulia/src/julia/miniBUDE.jl directory.
```shell
module load julia
julia --project=CUDA -e 'import Pkg; Pkg.instantiate()' 
julia --project=CUDA src/CUDA.jl -p 2 -w 32 --deck ../../../data/bm2 --device 2
```



