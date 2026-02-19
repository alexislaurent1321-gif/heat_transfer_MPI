#define FMT_HEADER_ONLY
#include "mpi.h"
#include "param.h"
#include "solve.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <fmt/core.h>

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

// Numarical calculation

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

// Solver

void updateT(const std::vector<double>& T, std::vector<double>& T_plus, int Nx_ghost, int Ny_ghost) {
    for(int i = 1; i < Nx_ghost - 1; ++i) {
        for(int j = 1; j < Ny_ghost - 1; ++j) {
           
            int index = i * Ny_ghost + j;
            T_plus[index] = T[index] + p.kappa * p.dt * (
                (T[(i+1)*Ny_ghost + j] - 2*T[index] + T[(i-1)*Ny_ghost + j]) / (p.dx*p.dx) +
                (T[index + 1] - 2*T[index] + T[index - 1]) / (p.dy*p.dy)
            );
        }
    }   
}

double solve(int N, int rank, int nprocs, std::shared_ptr<double> err_max, std::shared_ptr<double> Tmax) {
    // Setup Cartesian Topology
    int dims[2] = {0, 0};
    MPI_Dims_create(nprocs, 2, dims);
    int periods[2] = {0, 0};
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);

    int coords[2];
    MPI_Cart_coords(cart_comm, rank, 2, coords);
    int rank_up, rank_down, rank_left, rank_right;
    MPI_Cart_shift(cart_comm, 0, 1, &rank_up, &rank_down);   
    MPI_Cart_shift(cart_comm, 1, 1, &rank_left, &rank_right);

    // Define Local Dimensions
    int Nx_local = N / dims[0];
    int Ny_local = N / dims[1];
    int Nx_ghost = Nx_local + 2;
    int Ny_ghost = Ny_local + 2;

    std::vector<double> T(Nx_ghost * Ny_ghost, 0.0);
    std::vector<double> T_plus(Nx_ghost * Ny_ghost, 0.0);

    // Initial Condition (t=0)
    T_ex_fill(T, 0.0, Nx_ghost, Ny_ghost, coords);

    // Define Column Type for Non-Contiguous Memory
    MPI_Datatype column_type;
    MPI_Type_vector(Nx_local, 1, Ny_ghost, MPI_DOUBLE, &column_type);
    MPI_Type_commit(&column_type);

    double start_time = MPI_Wtime();

    for (int iter = 0; iter < p.Nt; ++iter) {
        // Exchange Rows (Up/Down)
        MPI_Sendrecv(&T[1*Ny_ghost + 1], Ny_local, MPI_DOUBLE, rank_up, 0,
                     &T[(Nx_ghost-1)*Ny_ghost + 1], Ny_local, MPI_DOUBLE, rank_down, 0, cart_comm, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&T[(Nx_ghost-2)*Ny_ghost + 1], Ny_local, MPI_DOUBLE, rank_down, 1,
                     &T[0*Ny_ghost + 1], Ny_local, MPI_DOUBLE, rank_up, 1, cart_comm, MPI_STATUS_IGNORE);

        // Exchange Columns (Left/Right)
        MPI_Sendrecv(&T[1*Ny_ghost + 1], 1, column_type, rank_left, 2,
                     &T[1*Ny_ghost + (Ny_ghost-1)], 1, column_type, rank_right, 2, cart_comm, MPI_STATUS_IGNORE);
        MPI_Sendrecv(&T[1*Ny_ghost + (Ny_ghost-2)], 1, column_type, rank_right, 3,
                     &T[1*Ny_ghost + 0], 1, column_type, rank_left, 3, cart_comm, MPI_STATUS_IGNORE);

        // Apply Global Boundary Conditions (Avoids "frozen" edges)
        apply_boundaries(T, iter * p.dt, Nx_ghost, Ny_ghost, coords, dims);

        // Compute Next Step
        updateT(T, T_plus, Nx_ghost, Ny_ghost);
        T.swap(T_plus);
    }
    
    double end_time = MPI_Wtime();

    // Data Export (One file per process for 2D reconstruction)
    std::ofstream file(fmt::format("T_data_{}_{}.txt", coords[0], coords[1]));
    for(int i = 1; i < Nx_ghost - 1; i++) {
        for(int j = 1; j < Ny_ghost - 1; j++) {
            file << T[i * Ny_ghost + j] << (j == Ny_ghost - 2 ? "" : " ");
        }
        file << "\n";
    }
    file.close();

    std::ofstream processFile("Process_datas.txt");
    processFile << dims[0] << " " << dims[1];
    processFile.close();

    // Global Error Reduction
    double l_err = error_T(T, Nx_ghost, Ny_ghost, coords);
    double l_Tmax = max_T(T, Nx_ghost, Ny_ghost);
    
    MPI_Reduce(&l_err, err_max.get(), 1, MPI_DOUBLE, MPI_MAX, 0, cart_comm);
    MPI_Reduce(&l_Tmax, Tmax.get(), 1, MPI_DOUBLE, MPI_MAX, 0, cart_comm);

    MPI_Type_free(&column_type);
    return end_time - start_time;
}