import numpy as np
import matplotlib.pyplot as plt

L = 8
H = 4 #of domain
Vinf = 1
rho = 1.225  # kg/m^3
nu = 1.5e-5
dx = 0.0008
Nx = int(L/dx)
Ny = 400


# physical grid
r = 1.012
a = Ny * np.log(r)
x = np.linspace(0, L, Nx)
y = np.array([H * ((np.exp(a*i/Ny))-1)/(np.exp(a)-1) for i in range(Ny)])
dy = np.diff(y)


# Computational grid (uniform)
eta = np.zeros((Ny))
ff = np.zeros((Ny))
for j in range(Ny):
    eta[j] = j/(Ny-1)
    ff[j] = (np.exp(a)-1)/np.exp(a*eta[j])/a/H

dxi = dx/L
deta = 1/(Ny-1)


# BCs
u = Vinf * np.ones((Nx, Ny))  # free stream
u[:, 0] = 0  # no-slip at the wall
v = np.zeros((Nx, Ny))  # initialize v to zero


for i in range(Nx - 1):
    # Tridiagonal system
    A = np.zeros((Ny, Ny))  # Coefficient matrix
    B = np.zeros(Ny)        # rhs

    for j in range(0, Ny-2):
        f = ff[j]
        aj = nu*f**2 *dxi*L/deta**2
        bj = 2 * nu * f**2 *L*dxi / deta**2   + u[i, j+1]
        cj = nu*f**2 *dxi*L/deta**2
        dj = - (u[i,j+2]-u[i,j])* (v[i,j+1]*f*L*dxi+nu*f**2 *L*dxi)/2/deta + u[i,j+1]**2

        # Append matrix elements
        A[j+1, j] = -aj
        A[j+1, j+1] = bj
        A[j+1, j+2] = -cj
        B[j+1] = dj

    # BCs
    A[0, 0] = 1
    B[0] = 0

    A[Ny-1, Ny-1] = Vinf
    B[Ny-1] = Vinf

    u[i+1,:] = np.linalg.solve(A, B)
    for j in range(0, Ny - 2):
        f = ff[j]
        v[i + 1, j + 1] = v[i + 1, j] - 0.5 * deta / L / f * (u[i + 1, j + 1] - u[i, j + 1] + u[i + 1, j] - u[i, j]) / dxi
    #BCs
    v[i, 0] = 0  # No-slip
    v[i, -1] = 0  # Free stream



# Boundary layer thickness
bl_thckns = np.zeros(Nx)
for t in range(0, Nx):
    for k in range(0, Ny):
        if u[t, k] > 0.99 * Vinf:
            bl_thckns[t] = y[k]
            break

# plt.plot(x, bl_thckns)
# plt.xlabel("x")
# plt.ylabel("Boundary Layer Thickness")
# plt.show()

d1 = np.zeros(Nx)
d2 = np.zeros(Nx)
dudytau = np.zeros(Nx)
tw = np.zeros(Nx)
cf = np.zeros(Nx)


for j in range(0, Nx):
    for k in range(0, Ny-1):
        d1[j] += (1 - u[j, k] / Vinf) * dy[k]
        d2[j] += (u[j, k] / Vinf) * (1 - u[j,k] / Vinf) * dy[k]
        dudytau[j] = (u[j, 1] - u[j, 0]) / dy[0]
        tw[j] = nu * rho * dudytau[j]
        cf[j] = tw[j] / (0.5 * rho * Vinf**2)


# Return function
def bon_implicit():
    return x, y, bl_thckns, d1, d2, cf, u, v
