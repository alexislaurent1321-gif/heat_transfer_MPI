#include "param.h"

/** \file    param.cpp
 *  \brief   Implementation of the functions to import parameters from a JSON file and calculate numerical parameters according to system size
 */

void Param::load(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return;

    json data;
    file >> data; 

    kappa   = data["physics"]["kappa"];
    t_final = data["physics"]["t_final"];
    Tmax    = data["physics"]["Tmax"];
    sigma   = data["physics"]["sigma"];
    
    sizes   = data["simulation"]["sizes"].get<std::vector<int>>();
}

// Calculation of numerical parameters according to system size
void Param::update(int size){
    dx = (xmax-xmin)/size;
    dy = (ymax-ymin)/size;
    dt = 0.5 * 1./(4 * kappa * (1/pow(dx,2) + 1/pow(dy,2)));        // dt is equal to half of the CFL condition
    Nt = static_cast<int>(t_final/dt);                              // number of iterations for finite differences
}