
#include <iostream>
#include <cmath>
#include <vector>



int main(){

    double high_energy;
    double low_energy;
    double p_energy;

    high_energy = 0.546031;
    low_energy = 0.545855;
    p_energy = 0.545872;

    double f = (high_energy - p_energy) / (high_energy - low_energy);

    std::cout << f << std::endl;

    return 0;

}