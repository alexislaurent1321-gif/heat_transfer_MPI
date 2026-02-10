#include <vector>
#include <memory>
#include <fstream>

#ifndef SOLVE_H
#define SOLVE_H

// Filling the T matrix with the exact solution
void T_ex(std::vector<double>& T, double t, int Nx_g, int Ny_g, int coords[2]);

// Exact solution of T evaluated at a point
double T_ex(double t, double x, double y); 

// Calculation of the infinite standard error between T and T_ex
double error_T(const std::vector<double>& T, int Nx_g, int Ny_g, int coords[2]);

// Calculation of the maximum value of T
double max_T(const std::vector<double>& T, int Nx_local, int Ny_local);

/*
Functions implementing the parallelized finite difference method
*/

// Update T using finite differences over a time iteration
void updateT(const std::vector<double>& T, std::vector<double>& T_plus, int Nx_ghost, int Ny_ghost, int coords[2]) ;   
  
// Resolution of T with MPI
double solve(int N, int rank, int nprocs, std::shared_ptr<double> err_max, std::shared_ptr<double> Tmax);

#endif