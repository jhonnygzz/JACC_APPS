#pragma once

#define MINIBUDE_VERSION ""
#define MINIBUDE_COMPILE_COMMANDS {  "/opt/nvidia/hpc_sdk/Linux_x86_64/24.9/compilers/bin/nvcc -forward-unknown-to-host-compiler -DCUDA -DUSE_PPWI=\"1\\\\,2\\\\,4\\\\,8\\\\,16\\\\,32\\\\,64\\\\,128\" -I<SRC>/cuda -I/home/xa2/Fall2024/JACC_APPS/miniBUDE/miniBUDECJulia/build_cudaH100/generated  -std=c++17 -forward-unknown-to-host-compiler -extended-lambda -use_fast_math -restrict -arch=sm_90   -DNDEBUG -O3 -march=native -ffast-math -std=c++17 -x cu -c <SRC>/main.cpp -o <OUT>/src/main.cpp.o" }
