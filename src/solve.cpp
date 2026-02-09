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
void updateT(const std::vector<double>& T, std::vector<double>& T_plus, int Nx_ghost, int Ny_ghost) {

    for(int i = 1; i < Nx_ghost - 1; ++i) {
        for(int j = 1; j < Ny_ghost - 1; ++j) {
            int idx = i * Ny_ghost + j;
            T_plus[idx] = T[idx] + p.kappa * p.dt * (
                (T[(i+1)*Ny_ghost + j] - 2*T[idx] + T[(i-1)*Ny_ghost + j]) / (p.dx*p.dx) +
                (T[idx + 1] - 2*T[idx] + T[idx - 1]) / (p.dy*p.dy)
            );
        }
    }   
}

// Function to execute the Jacobi method with MPI
double solve(int N, int rank, int nprocs, double *err_max, double *Tmax) {
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

    int Nx_local = N / dims[0];
    int Ny_local = N / dims[1];

    int Nx_ghost = Nx_local + 2;
    int Ny_ghost = Ny_local + 2;

    std::vector<double> T(Nx_ghost * Ny_ghost, 0.0);
    std::vector<double> T_plus(Nx_ghost * Ny_ghost, 0.0);

    // Creation of column type for MPI communication
    MPI_Datatype column_type;
    MPI_Type_vector(Nx_local, 1, Ny_ghost, MPI_DOUBLE, &column_type);
    MPI_Type_commit(&column_type);

    double start_time = MPI_Wtime();

    for (int iter = 0; iter < p.Nt; ++iter) {

        MPI_Sendrecv(&T[1*Ny_ghost + 1], Ny_local, MPI_DOUBLE, rank_up, 0,
                     &T[(Nx_ghost-1)*Ny_ghost + 1], Ny_local, MPI_DOUBLE, rank_down, 0, cart_comm, MPI_STATUS_IGNORE);
        
        MPI_Sendrecv(&T[(Nx_ghost-2)*Ny_ghost + 1], Ny_local, MPI_DOUBLE, rank_down, 1,
                     &T[0*Ny_ghost + 1], Ny_local, MPI_DOUBLE, rank_up, 1, cart_comm, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&T[1*Ny_ghost + 1], 1, column_type, rank_left, 2,
                     &T[1*Ny_ghost + (Ny_ghost-1)], 1, column_type, rank_right, 2, cart_comm, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&T[1*Ny_ghost + (Ny_ghost-2)], 1, column_type, rank_right, 3,
                     &T[1*Ny_ghost + 0], 1, column_type, rank_left, 3, cart_comm, MPI_STATUS_IGNORE);

        updateT(T, T_plus, Nx_ghost, Ny_ghost);
        T.swap(T_plus);
    }
    
    double end_time = MPI_Wtime();

    MPI_Type_free(&column_type);

    std::ofstream file(fmt::format("T_data_{}_{}.txt", coords[0], coords[1]));
    for(int i = 1; i < Nx_ghost - 1; i++) {
        for(int j = 1; j < Ny_ghost - 1; j++) {
            file << T[i * Ny_ghost + j] << (j == Ny_ghost - 2 ? "" : " ");
        }
        file << "\n";
    }
    file.close();

    *err_max = error_T(T, Nx_local, Ny_local, rank, nprocs);
    *Tmax = max_T(T, Nx_local, Ny_local, rank, nprocs);

    return end_time - start_time;
}

