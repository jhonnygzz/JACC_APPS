#include <iostream>
#include <hip/hip_runtime.h>

int main() {
    int deviceCount = 0;

    // Get the number of available GPUs
    hipError_t status = hipGetDeviceCount(&deviceCount);

    if (status != hipSuccess) {
        std::cerr << "Error: Unable to retrieve device count. HIP Error: " 
                  << hipGetErrorString(status) << std::endl;
        return 1;
    }

    if (deviceCount == 0) {
        std::cout << "No AMD GPUs found.\n";
    } else {
        std::cout << "Number of AMD GPUs found: " << deviceCount << std::endl;
    }

    return 0;
}