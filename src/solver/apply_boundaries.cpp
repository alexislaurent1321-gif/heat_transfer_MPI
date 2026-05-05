#include "mpi.h"
#include "param.h"
#include "solver.h"

/** \file    apply_boundaries.cpp
 *  \brief   Implementation of the function to apply Dirichlet boundary conditions to the global domain edges
 */

extern Param p;

// Apply Dirichlet boundary conditions to the global domain edges
void apply_boundaries(std::vector<double>& T, double t, int Nx_ghost, int Ny_ghost, int coords[2], int dims[2]) {
    int Nx_local = Nx_ghost - 2;
    int Ny_local = Ny_ghost - 2;

    // Global top border (i = 0)
    if (coords[0] == 0) {
 
        for (int j = 0; j < Ny_ghost; ++j) {
            double x = p.xmin - p.dx; 
            double y = p.ymin + (coords[1] * Ny_local + (j - 1)) * p.dy;
            
            T[0 * Ny_ghost + j] = T_ex(t, x, y);
        }
    }

    // Global bottom border (i = Nx_ghost - 1)
    if (coords[0] == dims[0] - 1) {
        
        for (int j = 0; j < Ny_ghost; ++j) {
            double x = p.xmax + p.dx;
            double y = p.ymin + (coords[1] * Ny_local + (j - 1)) * p.dy;
            
            T[(Nx_ghost - 1) * Ny_ghost + j] = T_ex(t, x, y);
        }
    }

    // Global left border (j = 0)
    if (coords[1] == 0) {
        
        for (int i = 0; i < Nx_ghost; ++i) {
            
            double x = p.xmin + (coords[0] * Nx_local + (i - 1)) * p.dx;
            double y = p.ymin - p.dy;
            
            T[i * Ny_ghost + 0] = T_ex(t, x, y);
        }
    }

    // Global right border (j = Ny_ghost - 1)
    if (coords[1] == dims[1] - 1) {

        for (int i = 0; i < Nx_ghost; ++i) {
            double x = p.xmin + (coords[0] * Nx_local + (i - 1)) * p.dx;
            double y = p.ymax + p.dy;

            T[i * Ny_ghost + (Ny_ghost - 1)] = T_ex(t, x, y);
        }
    }
}
