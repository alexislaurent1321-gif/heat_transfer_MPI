import numpy as np
import matplotlib.pyplot as plt

# Loading data

data1 = np.loadtxt("mpi_results_1.txt")
data2 = np.loadtxt("mpi_results_2.txt")
data4 = np.loadtxt("mpi_results_4.txt")


# Plot time

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


# Plot errors

error_data1 = np.loadtxt("error_1.txt")
error_data2 = np.loadtxt("error_2.txt")
error_data4 = np.loadtxt("error_4.txt")

error1 = error_data1[:,1]
error2 = error_data2[:,1]
error4 = error_data4[:,1]

plt.figure(figsize=(12,8))
plt.plot(np.log2(sizes), np.log2(error1), ':o', label=f'nprocs=1')
plt.plot(np.log2(sizes), np.log2(error2), '--', label=f'nprocs=2')
plt.plot(np.log2(sizes), np.log2(error4), ':+', label=f'nprocs=4')

plt.xlabel("log2(N)")
plt.ylabel("log2(Error)")

plt.legend()
plt.grid()
plt.savefig("graph/img/errors_t=5e4.jpg")
plt.show()

