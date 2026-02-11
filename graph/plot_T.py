import numpy as np
import matplotlib.pyplot as plt

def reconstruct_2d():
    full_matrix = []
    nprocs_x, nprocs_y = np.loadtxt("Process_datas.txt", dtype=int)

    for i in range(nprocs_x):
        row_blocks = []

        for j in range(nprocs_y):
            data = np.loadtxt(f"T_data_{i}_{j}.txt")
            row_blocks.append(data)

        full_matrix.append(np.hstack(row_blocks))
        
    return np.vstack(full_matrix)

# Plot
fig, ax = plt.subplots(figsize=(10,7))
T = reconstruct_2d()
im = ax.pcolormesh(T, cmap='rainbow', vmin=0, vmax=1350)

ax.set_xlabel("Nx")
ax.set_ylabel("Nz")

fig.colorbar(im, ax=ax)
plt.savefig("graph/img/T_t=5e4.jpg")

plt.show()
