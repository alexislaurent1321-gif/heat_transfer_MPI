#include <vector>
#include <cmath>

#ifndef PARAM_H
#define PARAM_H

/**
 * \brief Structure contenant les paramètres de la simulation
 */
struct Param{

    // Physical parameters

    // domain size
    double xmin=-0.5; 
    double xmax=0.5;
    double ymin=-0.5; 
    double ymax=0.5; 

    double kappa = 1e-6;    // constant diffusion coefficient
    double Tmax = 1350;     // amplitude
    double sigma = 0.3;     // standard deviation of the Gaussian solution
    double t_final = 1e4;   // time at which T is calculated


    // Numerical parameters

    std::vector<int> sizes = {64, 128, 256, 512}; // sizes of matrices to be tested
    double dx, dy, dt;
    int Nt;

    // Calculation of numerical parameters according to system size
    void update(int size){
        dx = (xmax-xmin)/size;
        dy = (ymax-ymin)/size;
        dt = 1./(4 * kappa * (1/pow(dx,2) + 1/pow(dy,2)));  // dt is equal to half the CFL condition
        Nt = (int)(t_final/dt);                             // number of iterations for finite differences
    }
};

#endif