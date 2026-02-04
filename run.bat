@echo off
:: Compilation of the C++ code
@REM mpic++ -std=c++20 -o test.exe main.cpp solve.cpp
g++ main.cpp solve.cpp -o test.exe -I"%MSMPI_INC%" -L"%MSMPI_LIB64%" -lmsmpi

:: Running with mpirun for 1, 2, and 4 processes
mpiexec -n 1 test.exe
mpiexec -n 2 test.exe
mpiexec -n 4 test.exe

:: Running Python scripts
python plot_T.py
python plot_mpi.py

:: Deletion of data files
del test.exe

del mpi_results_1.txt
del mpi_results_2.txt
del mpi_results_4.txt

del T_data_0.txt
del T_data_1.txt
del T_data_2.txt
del T_data_3.txt
del T_data.txt