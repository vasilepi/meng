# CFD 1st Assignment
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

L = 8
Vinf = 1
nu = 1.5e-5
Nx = 800
Ny = 300

dx = L/Nx
dy = 0.001

u = np.zeros((Ny, Nx))
v = np.zeros((Ny, Nx))

# BCs
u[0,:] = 0 # no-slip at wall
u[:,0] = Vinf # free stream
u[-1,:] = Vinf # free stream at top
v[0,:] = 0 # no-slip at wall
v[:,0] = 0 # free stream
v[-1,:] = 0 # free stream at top

bl_thckns = np.zeros(Nx)

for j in range(1, Nx):
    for i in range(1, Ny-1):
        dudy = (u[i+1,j-1] - u[i-1,j-1]) / (2 * dy)
        d2udy2 = (u[i+1,j-1] - 2 * u[i,j-1] + u[i-1, j-1]) / (dy ** 2)

        u[i, j] = u[i,j-1] + nu * d2udy2/u[i,j-1] * dx - v[i, j - 1] * dudy * dx / u[i,j-1]

    for z in range(Ny):
        if u[z, j] > 0.99 * Vinf:
            bl_thckns[j] = z * dy  # Store the y-height
            break

x = np.linspace(0, L, Nx)
plt.plot(x, bl_thckns, label="Boundary Layer Thickness")
plt.xlabel("x (m)")
plt.ylabel("Boundary Layer Thickness δ (m)")
plt.title("Boundary Layer Thickness over Flat Plate")
plt.legend()
plt.grid(True)
plt.show()

# # Visualize velocity profile as well
# plt.contourf(x, y, u, levels=50, cmap="viridis")
# plt.colorbar(label="Velocity (u)")
# plt.xlabel("x (m)")
# plt.ylabel("y (m)")
# plt.title("Boundary Layer Velocity Profile")
# plt.show()