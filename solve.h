#include <vector>
#include <fstream>

#ifndef SOLVE_H
#define SOLVE_H

using namespace std;

/**
 * \brief Filling the T matrix with the exact solution at a time t
 * 
 * \param T empty array for the temperature
 * \param t temperature calculated at time t
 * \param Nx spatial discretization on the x-axis
 * \param Ny spatial discretization on the y-axis
 * \param rank rank of the process
 * \param nprocs total number of processes
 */
void T_ex(vector<vector<double>>& T, double t, int Nx, int Ny, int rank, int nprocs);

/**
 * \brief Calculation of the exact T at a time t on a single point (useful for error calculation)
 * 
 * \param t temperature calculated at time t
 * \param x x coordinate of the point
 * \param y y coordinate of the point
 */
double T_ex(double t, double x, double y); 

/**
 * \brief Calculation of the error 
 * 
 * \param t temperature calculated at time t
 * \param Nx spatial discretization on the x-axis
 * \param Ny spatial discretization on the y-axis
 * \param rank rank of the process
 * \param nprocs total number of processes
 */
double error_T(vector<vector<double>> T, int Nx, int Ny, int rank, int nprocs);

/**
 * \brief Calculation of the maximal temperature in the domain at time t (useful for error calculation)
 * 
 * \param Nx spatial discretization on the x-axis
 * \param Ny spatial discretization on the y-axis
 * \param rank rank of the process
 * \param nprocs total number of processes
 */
double max_T(vector<vector<double>> T, int Nx, int Ny, int rank, int nprocs);



/*
Functions implementing the parallelized finite difference method
*/

/**
 * \brief Update T using finite differences over a time iteration
 * 
 * \param T last temperature
 * \param T_plus current temperature
 * \param Nx spatial discretization on the x-axis
 * \param Ny spatial discretization on the y-axis
 * \param rank rank of the process
 * \param nprocs total number of processes
 */
void updateT(const vector<vector<double>>& T, vector<vector<double>>& T_plus, int Nx, int Ny, int rank, int nprocs);   

/**
 * \brief Compute the new temperature by writing the temperature to the file corresponding to its process number
 * 
 * \param N number of points of the domain
 * \param rank rank of the process
 * \param nprocs total number of processes
 * \param err_max infinite standard error
 * \param Tmax maximale temperature on the domain (useful for error calculation)
 */
double solve(int N, int rank, int nprocs, double *err_max, double *Tmax);

#endif