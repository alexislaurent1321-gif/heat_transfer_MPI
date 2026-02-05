#define FMT_HEADER_ONLY
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>
#include "solve.h"
#include "param.h"
#include "merge.h"
#include <fmt/core.h>


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

    for (int N : p.sizes) {
        if(rank==0) cout << "N = " << N << endl;
        p.update(N); // Calculation of numerical parameters

        // Calculation of the infinite standard error for each thread
        double err_max_local = 0;
        double Tmax_local = 0;
        double parallel_time = solve(N, rank, nprocs, &err_max_local, &Tmax_local);

        // Maximum execution time for a thread
        double max_parallel_time;
        MPI_Reduce(&parallel_time, &max_parallel_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        
        // maximum error of T
        double err_max;
        MPI_Reduce(&err_max_local, &err_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        double Tmax;
        MPI_Reduce(&Tmax_local, &Tmax, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);


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

        vector<string> fileNames(nprocs);
        for(int rank=0; rank<nprocs; rank++){
            fileNames[rank] = fmt::format("T_data_{}.txt", rank);
        }
        mergeFiles(fileNames, "T_data.txt");    
    }
    
    MPI_Finalize();

    return 0;
}