# CFD 1st Assignment - explicit
import numpy as np
import matplotlib.pyplot as plt

L = 8
Vinf = 1
rho = 1.225 #kg/m^3
nu = 1.5e-5
Nx = 3000
Ny = 200

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

    # boundary layer thickness
    for z in range(Ny):
        if u[z, j] > 0.99 * Vinf:
            bl_thckns[j] = z * dy
            break

x = np.linspace(0, L, Nx)


d1 = np.zeros(Nx)
d2 = np.zeros(Nx)
dudytau = np.zeros(Nx)
tw = np.zeros(Nx)
cf = np.zeros(Nx)
for w in range(0, Nx):
    d1[w] = np.sum((1-u[:,w]/Vinf)*dy)
    d2[w] = np.sum((u[:, w] / Vinf) * (1 - u[:, w] / Vinf) * dy)
    dudytau[w] = (u[2,w]-u[1,w])/dy
    tw[w] = nu*rho*dudytau[w]
    cf[w] = tw[w]/(0.5*rho*Vinf**2)



def run_explicit():
    return x, bl_thckns, d1, d2, cf