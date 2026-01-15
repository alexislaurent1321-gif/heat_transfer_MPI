mpic++ -std=c++20 -o test.exe main.cpp solve.cpp

# Running with mpirun for 1, 2, and 4 processes
mpirun --use-hwthread-cpus -np 1 test.exe
mpirun --use-hwthread-cpus -np 2 test.exe
mpirun --use-hwthread-cpus -np 4 test.exe

# Running Python scripts
python3 plot_T.py
python3 plot_mpi.py

# Deletion of data files
rm test.exe

rm mpi_results_1.txt
rm mpi_results_2.txt
rm mpi_results_4.txt

rm T_data_0.txt
rm T_data_1.txt
rm T_data_2.txt
rm T_data_3.txt
rm T_data.txt


