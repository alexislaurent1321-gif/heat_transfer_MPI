mpic++ -std=c++20 -o test.exe main.cpp solve.cpp

# Exécution avec mpirun pour 1, 2 et 4 processus
mpirun --use-hwthread-cpus -np 1 test.exe
mpirun --use-hwthread-cpus -np 2 test.exe
mpirun --use-hwthread-cpus -np 4 test.exe

# Exécution des scripts Python
python3 plot_T.py
python3 plot_mpi.py

# Suppression des fichiers de données
rm test.exe

rm mpi_results_1.txt
rm mpi_results_2.txt
rm mpi_results_4.txt

rm T_data_0.txt
rm T_data_1.txt
rm T_data_2.txt
rm T_data_3.txt
rm T_data.txt


