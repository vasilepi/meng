# CFD 1st Assignment - implicit
import numpy as np
import matplotlib.pyplot as plt

L = 8
H = 0.1 #of domain
Vinf = 1
rho = 1.225  # kg/m^3
nu = 1.5e-5

dx = 0.0008
dy = 0.001
Nx = int(L/dx)
Ny = int(H/dy)

x = np.linspace(0, L, Nx)
y = np.linspace(0, Ny * dy, Ny)

# BCs
u = Vinf * np.ones((Ny, Nx))  # free stream
u[0, :] = 0  # no-slip at the wall
v = np.zeros((Ny, Nx))  # initialize v to zero


for j in range(Nx - 1):
    # Tridiagonal system
    A = np.zeros((Ny, Ny))  # Coefficient matrix
    B = np.zeros(Ny)        # rhs

    for i in range(1, Ny-1):
        aj = nu / dy**2
        bj = 2 * nu / dy**2 + u[i, j] / dx
        cj = nu / dy**2
        dj = -v[i, j] * (u[i+1, j] - u[i-1, j]) / (2 * dy) + u[i, j]**2 / dx

        # Append matrix elements
        A[i, i-1] = -aj
        A[i, i] = bj
        A[i, i+1] = -cj
        B[i] = dj

    # BCs
    A[0, 0] = 1
    B[0] = 0

    A[Ny-1, Ny-1] = Vinf
    B[Ny-1] = Vinf

    u[:,j+1] = np.linalg.solve(A, B)
    for i in range(1, Ny - 1):
        v[i, j + 1] = v[i-1, j - 1] - (dy / (2 * dx)) * ((u[i, j+1] - u[i, j]) + (u[i -1, j+1] - u[i - 1, j]))

    #BCs
    v[0, j] = 0  # No-slip
    v[-1, j] = 0  # Free stream





# Boundary layer thickness
bl_thckns = np.zeros(Nx)
for j in range(Nx):
    bl_thckns[j] = y[np.where(u[:, j] > 0.99 * Vinf)[0][0]]

d1 = np.zeros(Nx)
d2 = np.zeros(Nx)
dudytau = np.zeros(Nx)
tw = np.zeros(Nx)
cf = np.zeros(Nx)

for j in range(0, Nx):
    d1[j] = np.sum((1 - u[:, j] / Vinf) * dy)
    d2[j] = np.sum((u[:, j] / Vinf) * (1 - u[:, j] / Vinf) * dy)
    dudytau[j] = (u[2, j] - u[1, j]) / dy  # Gradient near the wall
    tw[j] = nu * rho * dudytau[j]
    cf[j] = tw[j] / (0.5 * rho * Vinf**2)

# Return function
def run_implicit():
    return x, y, bl_thckns, d1, d2, cf, u, v
