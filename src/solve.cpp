#define FMT_HEADER_ONLY
#include "mpi.h"
#include "param.h"
#include "solve.h"
#include <cmath>
#include <iostream>
#include <fmt/core.h>

extern Param p; // Importing settings

// Functions related to the exact temperature

/* Exact solution for temperature
Filling T at time t
Block division on x
*/
void T_ex(std::vector<double>& T, double t, int Nx_local, int Ny, int rank, int nprocs){
    
    for(int i=0; i < Nx_local; ++i){
        double x = p.xmin + (rank * (Nx_local - 2) + (i - 1)) * p.dx;
        
        for(int j=0; j < Ny; ++j){
            double y = p.ymin + j * p.dy;
            T[i * Ny + j] = T_ex(t, x, y);
        }
    }
}

// Exact solution of T evaluated at a point
double T_ex(double t, double x, double y){
    double denom = 2 * pow(p.sigma, 2) + 4 * t * p.kappa;
    double ampl = p.Tmax / (1 + 2 * t * p.kappa / pow(p.sigma, 2));
    return ampl * exp((-pow(x, 2) - pow(y, 2)) / denom);
}

// Calculation of the infinite standard error between T and T_ex
double error_T(const std::vector<double>& T, int Nx, int Ny, int rank, int nprocs){
    
    double err = 0;
    
    for(int i=1; i < Nx-1; i++){
        double x = p.xmin + (rank * (Nx - 2) + (i - 1)) * p.dx;
        
        for(int j=0; j < Ny; j++){
            double y = p.ymin + j * p.dy;
            double err_ij = std::abs(T[i * Ny + j] - T_ex(p.t_final, x, y));
            if(err < err_ij) err = err_ij;
        }
    }
    
    return err;
}

// Calculation of the infinite standard error between T and T_ex
double max_T(std::vector<double>& T, int Nx, int Ny, int rank, int nprocs){

    double Tmax = 0;
    double T_ij;    // error evaluated at a point
    double x, y;    // definition of x and y based on indices i and j

    for(int i=1; i<Nx-1; i++){
        x = p.xmin + (rank*Nx+i)*(p.xmax-p.xmin)/(Nx*nprocs);

        for(int j=0; j<Ny; j++){
            y = p.ymin + j*(p.ymax-p.ymin)/Ny;
            T_ij = T_ex(p.t_final,x,y);
            if(fabs(Tmax)<fabs(T_ij)) Tmax=T_ij;
        }
    }

    return Tmax;
}




/*
Numerical resolution functions
*/

// Update T using finite differences
void updateT(const std::vector<double>& T, std::vector<double>& T_plus, int Nx, int Ny) {
    int index;
    for(int i=1; i < Nx-1; ++i) {
        
        for(int j=1; j < Ny-1; ++j) {
            index = i * Ny + j;
            // Correction des parenthèses d'indexation
            T_plus[index] = T[index] + p.kappa * p.dt * (
                (T[(i+1)*Ny + j] - 2*T[index] + T[(i-1)*Ny + j]) / (p.dx*p.dx) +
                (T[index + 1] - 2*T[index] + T[index - 1]) / (p.dy*p.dy)
            );
        }
    }   
}


// Function to execute the Jacobi method with MPI
double solve(int N, int rank, int nprocs, double *err_max, double *Tmax) {
    // N_local inclut 2 lignes fantômes (sauf bords extrêmes, mais simplifié ici)
    int Nx_local = N / nprocs + 2;
    int Ny = N;

    std::vector<double> T(Nx_local * Ny, 0.0);
    std::vector<double> T_plus(Nx_local * Ny, 0.0);

    T_ex(T, 0, Nx_local, Ny, rank, nprocs);
    T_plus = T;

    double start_time = MPI_Wtime();

    for (int iter = 0; iter < p.Nt; ++iter) {
        // Échanges des Ghost Cells (Bloquants mais sûrs pour débuter)
        if (rank > 0) {
            MPI_Sendrecv(&T[Ny], Ny, MPI_DOUBLE, rank-1, 0, 
                         &T[0], Ny, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < nprocs-1) {
            MPI_Sendrecv(&T[(Nx_local-2)*Ny], Ny, MPI_DOUBLE, rank+1, 0, 
                         &T[(Nx_local-1)*Ny], Ny, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        updateT(T, T_plus, Nx_local, Ny);
        T.swap(T_plus);
    }
    
    double end_time = MPI_Wtime();

    // Sauvegarde propre
    std::ofstream file(fmt::format("T_data_{}.txt", rank));
    for(int i=1; i < Nx_local-1; i++){
        for(int j=0; j < Ny; j++){
            file << T[i * Ny + j] << (j == Ny-1 ? "" : " ");
        }
        file << "\n";
    }
    file.close();

    *err_max = error_T(T, Nx_local, Ny, rank, nprocs);
    *Tmax = max_T(T, Nx_local, Ny, rank, nprocs);

    return end_time - start_time;
}