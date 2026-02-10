#include <mpi.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>
#include "solve.h"
#include "param.h"


Param p; // Importing settings
using namespace std;

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    // Files for saving results
    ofstream results;  // file containing performance results
    if (rank == 0) {
        cout << "nprocs = " << nprocs << endl;
        if(nprocs==1) results.open("mpi_results_1.txt");
        if(nprocs==2) results.open("mpi_results_2.txt");
        if(nprocs==4) results.open("mpi_results_4.txt");
    }

    p.load("parameters.json"); // Load settings from config file

    for (int N : p.sizes) {
        if(rank==0) cout << "N = " << N << endl;
        p.update(N); // Calculation of numerical parameters

        // Calculation of the infinite standard error for each thread
        std::shared_ptr<double> err_max_local = std::make_shared<double>(0.);
        std::shared_ptr<double> Tmax_local = std::make_shared<double>(0.);
        double parallel_time = solve(N, rank, nprocs, err_max_local, Tmax_local);

        // Maximum execution time for a thread
        double max_parallel_time;
        MPI_Reduce(&parallel_time, &max_parallel_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        
        // maximum error of T
        double err_max = *err_max_local;
        double Tmax = *Tmax_local;

        // relative max error (%)
        double err_max_rel = err_max / Tmax * 100.; 

        if (rank == 0) {
            results << N << "  " << max_parallel_time << endl;
            cout << "Erreur relative en norme infinie = " << err_max_rel << "%" << endl;
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) {
        results.close();
        cout << "Résultats sauvegardés dans 'mpi_results.txt'.\n";
    }
    
    MPI_Finalize();

    return 0;
}