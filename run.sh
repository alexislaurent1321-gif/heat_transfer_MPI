mkdir -p build
cd build
cmake -DCMAKE_CXX_COMPILER=mpicxx ..
make
cd ..

mpirun --use-hwthread-cpus -np 1 ./build/heat_solver
mpirun --use-hwthread-cpus -np 2 ./build/heat_solver
mpirun --use-hwthread-cpus -np 4 ./build/heat_solver


if [ -d ".venv/" ]; then
    . .venv/bin/activate
fi

python3 graph/plot_T.py
python3 graph/plot_mpi.py

find . -maxdepth 1 -type f -name "*.txt" ! -name "CMakeLists.txt" -delete
