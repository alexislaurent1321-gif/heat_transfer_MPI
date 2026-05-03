#include "mpi.h"
#include "param.h"
#include "solver.h"

extern Param p; 

/*
Exact solution
*/

// Exact solution evaluated at a point (x,y) and time t
double T_ex(double t, double x, double y) {
    double denom = 2 * pow(p.sigma, 2) + 4 * t * p.kappa;
    double ampl = p.Tmax / (1 + 2 * t * p.kappa / pow(p.sigma, 2));
    return ampl * exp((-pow(x, 2) - pow(y, 2)) / denom);
}


// Fill the local block (including ghost cells) with the exact solution
void T_ex_fill(std::vector<double>& T, double t, int Nx_ghost, int Ny_ghost, int coords[2]) {
    // Use of ghost cells 
    int Nx_local = Nx_ghost - 2;    
    int Ny_local = Ny_ghost - 2;

    for(int i = 0; i < Nx_ghost; ++i) {
        double x = p.xmin + (coords[0] * Nx_local + (i - 1)) * p.dx;

        for(int j = 0; j < Ny_ghost; ++j) {
            double y = p.ymin + (coords[1] * Ny_local + (j - 1)) * p.dy;
            
            T[i * Ny_ghost + j] = T_ex(t, x, y);
        }
    }
}
