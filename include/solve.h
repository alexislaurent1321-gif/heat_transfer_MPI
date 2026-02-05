#include <vector>
#include <fstream>

#ifndef SOLVE_H
#define SOLVE_H

// Filling the T matrix with the exact solution
void T_ex(std::vector<double>& T, double t, int Nx, int Ny, int rank, int nprocs); 

// Exact solution of T evaluated at a point
double T_ex(double t, double x, double y); 

// Calculation of the infinite standard error between T and T_ex
double error_T(const std::vector<double>& T, int Nx, int Ny, int rank, int nprocs);

// Calculation of the maximum value of T
double max_T(const std::vector<double>& T, int Nx, int Ny, int rank, int nprocs);

/*
Functions implementing the parallelized finite difference method
*/

// Update T using finite differences over a time iteration
void updateT(const std::vector<double>& T, std::vector<double>& T_plus, int Nx, int Ny);   
  
// Resolution of T with MPI
double solve(int N, int rank, int nprocs, double *err_max, double *Tmax);

#endif