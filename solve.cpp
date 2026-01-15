#define FMT_HEADER_ONLY
#include "mpi.h"
#include "Param.h"
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
void T_ex(vector<vector<double>>& T, double t, int Nx, int Ny, int rank, int nprocs){
    double x, y;
    for(int i=0; i<Nx; ++i){
        x = p.xmin + (rank*Nx+i)*(p.xmax-p.xmin)/(Nx*nprocs);
        for(int j=0; j<Ny; ++j){
            y = p.ymin + j*(p.ymax-p.ymin)/Ny;
            T[i][j] = p.Tmax/(1+2*t*p.kappa/pow(p.sigma,2))*exp((-pow(x,2)-pow(y,2))/(2*pow(p.sigma,2)+4*t*p.kappa));
        }
    }
}

// Exact solution of T evaluated at a point
double T_ex(double t, double x, double y){
    return p.Tmax/(1+2*t*p.kappa/pow(p.sigma,2))*exp((-pow(x,2)-pow(y,2))/(2*pow(p.sigma,2)+4*t*p.kappa));
}

// Calculation of the infinite standard error between T and T_ex
double error_T(vector<vector<double>> T, int Nx, int Ny, int rank, int nprocs){
    double err = 0;     // max arror
    double err_ij;      // error evaluated at a point
    double x, y;        // definition of x and y based on indices i and j
    for(int i=1; i<Nx-1; i++){
        x = p.xmin + (rank*Nx+i)*(p.xmax-p.xmin)/(Nx*nprocs);
        for(int j=0; j<Ny; j++){
            y = p.ymin + j*(p.ymax-p.ymin)/Ny;
            err_ij = fabs(T[i][j]-T_ex(p.t_final,x,y));
            if(err<err_ij) err=err_ij;
        }
    }

    return err;
}

// Calculation of the infinite standard error between T and T_ex
double max_T(vector<vector<double>> T, int Nx, int Ny, int rank, int nprocs){
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
void updateT(const vector<vector<double>>& T, vector<vector<double>>& T_plus, int Nx, int Ny, int rank, int nprocs) {
    for(int i=1; i<Nx-1; ++i) {
        for(int j=1; j<Ny-1; ++j) {
            T_plus[i][j] = T[i][j] + p.kappa * p.dt * (
                (T[i+1][j] - 2*T[i][j] + T[i-1][j]) / (p.dx*p.dx) +
                (T[i][j+1] - 2*T[i][j] + T[i][j-1]) / (p.dy*p.dy));
        }
    }   
}


// Function to execute the Jacobi method with MPI
double solve(int N, int rank, int nprocs, double *err_max, double *Tmax) {

    int N_local = 0;
    if(nprocs==1) 
        N_local = N;
    else
        N_local = N / nprocs + 2; // Number of local lines (includes ghost lines)

    // Creation and filling of block T with the solution at initial time
    vector<vector<double>> T(N_local,vector<double>(N,0));
    vector<vector<double>> T_plus(N_local,vector<double>(N,0));
    T_ex(T,0,N_local,N,rank,nprocs);
    T_ex(T_plus,0,N_local,N,rank,nprocs);

    double start_time = MPI_Wtime();

    for (int iter = 0; iter <= p.Nt; ++iter) {
        // Exchange border values with neighboring threads
        if (rank > 0) {
            MPI_Send(&T[1][0], N, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
            MPI_Recv(&T[0][0], N, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < nprocs-1) {
            MPI_Recv(&T[N_local-1][0], N, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&T[N_local-2][0], N, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // Calculation of the update of T
        updateT(T, T_plus, N_local, N, rank, nprocs);
        T.swap(T_plus);
    }

    // Printing the table T
    std::vector<ofstream> T_file(nprocs);
    T_file[rank].open(fmt::format("T_data_{}.txt", rank));
    for(int i=1; i<N_local-1; i++){
        T_file[rank] << "\n";
        for(int j=0; j<N; j++){
            T_file[rank]<< T[i][j] << " ";
        }
    }
    T_file[rank].close();

    
    // Calculation of the inf error
    *err_max = error_T(T_plus, N_local, N, rank, nprocs);
    *Tmax = max_T(T_plus, N_local, N, rank, nprocs);

    double end_time = MPI_Wtime();
    return end_time - start_time;
}