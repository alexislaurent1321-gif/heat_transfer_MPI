@echo off

mpic++ -std=c++20 -o test.exe main.cpp solve.cpp

:: Exécution avec mpirun pour 1, 2 et 4 processus
mpiexec --use-hwthread-cpus -n 1 test.exe
mpiexec --use-hwthread-cpus -n 2 test.exe
mpiexec --use-hwthread-cpus -n 4 test.exe

:: Exécution des scripts Python
python plot_T.py
python plot_mpi.py

:: Suppression des fichiers de données
del test.exe

del mpi_results_1.txt
del mpi_results_2.txt
del mpi_results_4.txt

del T_data_0.txt
del T_data_1.txt
del T_data_2.txt
del T_data_3.txt
del T_data.txt
