mkdir -p build
cd build
cmake -DCMAKE_CXX_COMPILER=mpicxx ..
make
cd ..

mpirun --use-hwthread-cpus -np 1 ./build/heat_solver
mpirun --use-hwthread-cpus -np 2 ./build/heat_solver
mpirun --use-hwthread-cpus -np 4 ./build/heat_solver

python3 graph/plot_T.py
python3 graph/plot_mpi.py

rm mpi_results_*.txt 
rm T_data_*.txt T_data.txt