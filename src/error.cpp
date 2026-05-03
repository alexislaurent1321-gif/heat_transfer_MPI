#include "mpi.h"
#include "param.h"
#include "solve.h"

extern Param p; 

// Error

double error_T(const std::vector<double>& T, int Nx_ghost, int Ny_ghost, int coords[2]) {
    double err_max_local = 0;
    int Nx_local = Nx_ghost - 2;
    int Ny_local = Ny_ghost - 2;

    for(int i = 1; i < Nx_ghost - 1; ++i) {
        double x = p.xmin + (coords[0] * Nx_local + (i - 1)) * p.dx;

        for(int j = 1; j < Ny_ghost - 1; ++j) {
            double y = p.ymin + (coords[1] * Ny_local + (j - 1)) * p.dy;
            double diff = std::abs(T[i * Ny_ghost + j] - T_ex(p.t_final, x, y));

            if (diff > err_max_local) 
                err_max_local = diff;
        }
    }
    return err_max_local;
}

double max_T(const std::vector<double>& T, int Nx_ghost, int Ny_ghost) {
    double Tmax_local = 0;

    for(int i = 1; i < Nx_ghost - 1; ++i) {
        
        for(int j = 1; j < Ny_ghost - 1; ++j) {
            double val = std::abs(T[i * Ny_ghost + j]);
            
            if(val > Tmax_local) 
                Tmax_local = val;
        }
    }
    return Tmax_local;
}