#include <vector>
#include <cmath>
#include <nlohmann/json.hpp>
#include <fstream>

#ifndef PARAM_H
#define PARAM_H

using json = nlohmann::json;

/**
 * \brief Structure containing the simulation parameters
 */
struct Param{

    // Physical parameters

    // domain size
    double xmin=-0.5; 
    double xmax=0.5;
    double ymin=-0.5; 
    double ymax=0.5; 

    double kappa;       // constant diffusion coefficient
    double Tmax;        // amplitude
    double sigma;       // standard deviation of the Gaussian solution
    double t_final;     // time at which T is calculated


    // Numerical parameters

    std::vector<int> sizes; // sizes of matrices to be tested
    double dx, dy, dt;
    int Nt;

    void load(const std::string& filename) {
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
    void update(int size){
        dx = (xmax-xmin)/size;
        dy = (ymax-ymin)/size;
        dt = 0.5 * 1./(4 * kappa * (1/pow(dx,2) + 1/pow(dy,2)));  // dt is equal to half the CFL condition
        Nt = static_cast<int>(t_final/dt);                             // number of iterations for finite differences
    }
};

#endif