import numpy as np
import matplotlib.pyplot as plt


L = 8
H = 2     # of domain
Vinf = 1
rho = 1.225  # kg/m^3
nu = 1.5e-5
Ny = 300
dx = 0.0008
Nx = int(L / dx)


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
    eta[j] = a*j/Ny
    ff[j] = (np.exp(a)-1)/np.exp(eta[j])/H

dxi = dx/L
deta = a/Ny
u = np.zeros((Nx, Ny))
v = np.zeros((Nx, Ny))

# BCs
u[:,0] = 0 # no-slip at wall
u[0,:] = Vinf # free stream
u[:,-1] = Vinf # free stream at top

bl_thckns = np.zeros(Nx)

for i in range(0, Nx-1):
    for j in range(0, Ny-2):
        dudeta = 0.5*(u[i,j+2] - u[i,j])/deta
        d2udeta2 = (u[i,j+2]-2*u[i,j+1]+u[i,j])/deta**2
        f = ff[j]
        A = (nu*f**2 * L*dxi + v[i,j+1]*f*L*dxi)/u[i,j+1]
        u[i+1,j+1] = u[i,j+1] + nu*f**2 * d2udeta2 * L*dxi/u[i,j+1] - A*dudeta
        v[i+1,j+1] = v[i+1,j] - 0.5*deta/L/f *(u[i+1,j+1] - u[i,j+1] + u[i+1,j] - u[i,j]) /dxi

for t in range(0,Nx):
    for k in range(0,Ny):
        if u[t,k] > 0.99*Vinf:
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
def bon_explicit():
    return x, y, bl_thckns, d1, d2, cf, u, v
