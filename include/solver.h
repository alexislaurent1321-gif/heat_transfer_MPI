#include <vector>
#include <memory>
#include <fstream>

#ifndef SOLVE_H
#define SOLVE_H

/** \file    solver.h
 *  \brief   Header file for the solver functions
 */



/** \brief   Exact solution of T evaluated at a point
 *  \param t  Time
 *  \param x  X-coordinate
 *  \param y  Y-coordinate
 *  \return   Exact temperature at the given point and time
 */
double T_ex(double t, double x, double y); 


/** \brief   Fill the local block (including ghost cells) with the exact solution
 *  \param T    Vector to be filled
 *  \param t    Time
 *  \param Nx_g Number of cells in the x-direction (including ghosts)
 *  \param Ny_g Number of cells in the y-direction (including ghosts)
 *  \param coords Coordinates of the process in the Cartesian topology
 */
void T_ex_fill(std::vector<double>& T, double t, int Nx_g, int Ny_g, int coords[2]);


/** \brief   Calculate the infinite standard error between T and T_ex
 *  \param T      Vector containing the numerical solution
 *  \param Nx_g   Number of cells in the x-direction (including ghosts)
 *  \param Ny_g   Number of cells in the y-direction (including ghosts)
 *  \param coords Coordinates of the process in the Cartesian topology
 *  \return       Infinite standard error
 */
double error_T(const std::vector<double>& T, int Nx_g, int Ny_g, int coords[2]);


/** \brief   Calculate the maximum value of T
 *  \param T         Vector containing the numerical solution
 *  \param Nx_local  Number of cells in the x-direction (excluding ghosts)
 *  \param Ny_local  Number of cells in the y-direction (excluding ghosts)
 *  \return          Maximum value of T
 */
double max_T(const std::vector<double>& T, int Nx_local, int Ny_local);



/*
Functions implementing the parallelized finite difference method
*/


/** \brief   Apply Dirichlet boundary conditions to the global domain edges
 *  \param T      Vector containing the temperature
 *  \param t      Time
 *  \param Nx_g   Number of cells in the x-direction (including ghosts)
 *  \param Ny_g   Number of cells in the y-direction (including ghosts)
 *  \param coords Coordinates of the process in the Cartesian topology
 *  \param dims   Dimensions of the Cartesian topology
 */
void apply_boundaries(std::vector<double>& T, double t, int Nx_ghost, int Ny_ghost, int coords[2], int dims[2]);


/** \brief   Update T using finite differences over a time iteration
 *  \param T      Vector containing the current temperature
 *  \param T_plus Vector containing the updated temperature
 *  \param Nx_g   Number of cells in the x-direction (including ghosts)
 *  \param Ny_g   Number of cells in the y-direction (including ghosts)
 *  \param coords Coordinates of the process in the Cartesian topology
 */
void updateT(const std::vector<double>& T, std::vector<double>& T_plus, int Nx_ghost, int Ny_ghost, int coords[2]) ;   
  

/** \brief   Resolution of T with MPI
 *  \param N      Total number of cells in each direction
 *  \param rank   Rank of the current process
 *  \param nprocs Number of processes
 *  \param err_max Shared pointer to store the maximum error
 *  \param Tmax    Shared pointer to store the maximum temperature
 *  \return       Final temperature value
 */
double solve(int N, int rank, int nprocs, std::shared_ptr<double> err_max, std::shared_ptr<double> Tmax);

#endif