#include <vector>
#include <fstream>

#ifndef SOLVE_H
#define SOLVE_H

using namespace std;

void T_ex(vector<vector<double>>& T, double t, int Nx, int Ny, int rank, int nprocs); // Filling the T matrix with the exact solution

double T_ex(double t, double x, double y); // Exact solution of T evaluated at a point

// Calculation of the infinite standard error between T and T_ex
double error_T(vector<vector<double>> T, int Nx, int Ny, int rank, int nprocs);

double max_T(vector<vector<double>> T, int Nx, int Ny, int rank, int nprocs);



/*
Functions implementing the parallelized finite difference method
*/

// Update T using finite differences over a time iteration
void updateT(const vector<vector<double>>& T, vector<vector<double>>& T_plus, int Nx, int Ny, int rank, int nprocs);   

// Resolution of T with MPI
double solve(int N, int rank, int nprocs, double *err_max, double *Tmax);

#endif