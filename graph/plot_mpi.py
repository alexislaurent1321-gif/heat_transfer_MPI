import numpy as np
import matplotlib.pyplot as plt

# Loading data

data1 = np.loadtxt("mpi_results_1.txt")
data2 = np.loadtxt("mpi_results_2.txt")
data4 = np.loadtxt("mpi_results_4.txt")


# Plot

sizes = data1[:,0] 
time1 = data1[:,1]
time2 = data2[:,1]
time4 = data4[:,1]

plt.figure(figsize=(12,8))
plt.plot(sizes, time1 / time1, marker='o', label=f'nprocs=1')
plt.plot(sizes, time1 / time2, marker='o', label=f'nprocs=2')
plt.plot(sizes, time1 / time4, marker='o', label=f'nprocs=4')

plt.xlabel("N")
plt.ylabel("speedup")
plt.legend()
plt.grid()
plt.savefig("graph/img/performances_t=5e4.jpg")
plt.show()