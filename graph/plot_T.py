import numpy as np
import matplotlib.pyplot as plt

T = np.loadtxt("T_data.txt")

# Plot
fig, ax = plt.subplots(figsize=(10,7))
im = ax.pcolormesh(T, cmap='rainbow', vmin=0, vmax=1350)

ax.set_xlabel("Nx")
ax.set_ylabel("Nz")

fig.colorbar(im, ax=ax)
plt.savefig("T_t=5e4.jpg")

plt.show()
