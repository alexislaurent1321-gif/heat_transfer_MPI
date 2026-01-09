#include <vector>
#include <fstream>

#ifndef SOLVE_H
#define SOLVE_H

using namespace std;

void T_ex(vector<vector<double>>& T, double t, int Nx, int Ny, int rank, int nprocs); // Remplissage de la matrice T par la solution exacte

double T_ex(double t, double x, double y); // Solution exacte de T évaluée en un point

// Calcul de l'erreur en norme infinie entre T et T_ex
double error_T(vector<vector<double>> T, int Nx, int Ny, int rank, int nprocs);

double max_T(vector<vector<double>> T, int Nx, int Ny, int rank, int nprocs);



/*
Fonctions mettant en place la méthode des différences finies parallélisée
*/

// Mise à jour de T par différences finies sur une itération temporelle
void updateT(const vector<vector<double>>& T, vector<vector<double>>& T_plus, int Nx, int Ny, int rank, int nprocs);   

// Résolution de T avec MPI
double solve(int N, int rank, int nprocs, double *err_max, double *Tmax);

#endif