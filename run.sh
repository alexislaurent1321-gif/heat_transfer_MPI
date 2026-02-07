mkdir -p build
cd build
cmake -DCMAKE_CXX_COMPILER=mpicxx ..
make
cd ..

mpirun --use-hwthread-cpus -np 1 ./build/heat_solver
mpirun --use-hwthread-cpus -np 2 ./build/heat_solver
mpirun --use-hwthread-cpus -np 4 ./build/heat_solver
mpirun --use-hwthread-cpus -np 8 ./build/heat_solver

if [ -d ".venv/" ]; then
    . .venv/bin/activate
fi

python3 graph/plot_T.py
python3 graph/plot_mpi.py

rm mpi_results*.txt 
rm T_data*.txt 