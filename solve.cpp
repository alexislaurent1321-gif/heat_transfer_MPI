#define FMT_HEADER_ONLY
#include "mpi.h"
#include "param.h"
#include "solve.h"
#include <cmath>
#include <iostream>
#include <fmt/core.h>

extern Param p; // Importation des paramètres

// Fonctions concernant la température exacte

/* Solution exacte pour la température
Remplissage de T au temps t
Division par bloc sur x
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

// Solution exacte de T évaluée en un point
double T_ex(double t, double x, double y){
    return p.Tmax/(1+2*t*p.kappa/pow(p.sigma,2))*exp((-pow(x,2)-pow(y,2))/(2*pow(p.sigma,2)+4*t*p.kappa));
}

// Calcul de l'erreur en norme infinie entre T et T_ex
double error_T(vector<vector<double>> T, int Nx, int Ny, int rank, int nprocs){
    double err = 0; // erreur max
    double err_ij; // erreur évaluée en un point
    double x, y; // définition de x et y en fonction des indices i et j
    for(int i=1; i<Nx-1; i++){
        x = p.xmin + (rank*Nx+i)*(p.xmax-p.xmin)/(Nx*nprocs);
        for(int j=0; j<Ny; j++){
            y = p.ymin + j*(p.ymax-p.ymin)/Ny;
            err_ij = fabs(T[i][j]-T_ex(p.t_final,x,y));
            if(err<err_ij) err=err_ij;
        }
    }
    // std::cout << nprocs << std::endl;
    return err;
}

// Calcul de l'erreur en norme infinie entre T et T_ex
double max_T(vector<vector<double>> T, int Nx, int Ny, int rank, int nprocs){
    double Tmax = 0;
    double T_ij; // erreur évaluée en un point
    double x, y; // définition de x et y en fonction des indices i et j
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
Fonctions de résolution numérique
*/

// Mise à jour de T par différences finies
void updateT(const vector<vector<double>>& T, vector<vector<double>>& T_plus, int Nx, int Ny, int rank, int nprocs) {
    for(int i=1; i<Nx-1; ++i) {
        for(int j=1; j<Ny-1; ++j) {
            T_plus[i][j] = T[i][j] + p.kappa * p.dt * (
                (T[i+1][j] - 2*T[i][j] + T[i-1][j]) / (p.dx*p.dx) +
                (T[i][j+1] - 2*T[i][j] + T[i][j-1]) / (p.dy*p.dy));
        }
    }   
}


// Fonction pour exécuter la méthode Jacobi avec MPI
double solve(int N, int rank, int nprocs, double *err_max, double *Tmax) {
    // cout << p.dx << endl;
    int N_local = 0;
    if(nprocs==1) 
        N_local = N;
    else
        N_local = N / nprocs + 2; // Nombre de lignes locales (inclut lignes fantômes)

    // Création et remplissage du bloc T avec la solution au temps initial
    vector<vector<double>> T(N_local,vector<double>(N,0));
    vector<vector<double>> T_plus(N_local,vector<double>(N,0));
    T_ex(T,0,N_local,N,rank,nprocs);
    T_ex(T_plus,0,N_local,N,rank,nprocs);

    double start_time = MPI_Wtime();

    for (int iter = 0; iter <= p.Nt; ++iter) {
        // Échanger les valeurs des frontières avec les threads voisins
        if (rank > 0) {
            MPI_Send(&T[1][0], N, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
            MPI_Recv(&T[0][0], N, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < nprocs-1) {
            MPI_Recv(&T[N_local-1][0], N, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(&T[N_local-2][0], N, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
        }

        MPI_Barrier(MPI_COMM_WORLD);

        // Calcul de la mise à jour de T
        updateT(T, T_plus, N_local, N, rank, nprocs);
        T.swap(T_plus);
    }

    // Impression du tableau T
    std::vector<ofstream> T_file(nprocs);
    T_file[rank].open(fmt::format("T_data_{}.txt", rank));
    for(int i=1; i<N_local-1; i++){
        T_file[rank] << "\n";
        for(int j=0; j<N; j++){
            T_file[rank]<< T[i][j] << " ";
        }
    }
    T_file[rank].close();

    
    // Calcul de l'erreur inf
    *err_max = error_T(T_plus, N_local, N, rank, nprocs);
    *Tmax = max_T(T_plus, N_local, N, rank, nprocs);

    double end_time = MPI_Wtime();
    return end_time - start_time;
}