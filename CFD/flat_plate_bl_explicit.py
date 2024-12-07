# CFD 1st Assignment - explicit
import numpy as np
import matplotlib.pyplot as plt

L = 8
H = 0.1 #of domain
Vinf = 1
rho = 1.225 #kg/m^3
nu = 1.5e-5

dx = 0.0001
dy = 0.001
Nx = int(L/dx)
Ny = int(H/dy)

u = np.zeros((Nx, Ny))
v = np.zeros((Nx, Ny))

# BCs
u[:,0] = 0 # no-slip at wall
u[0,:] = Vinf # free stream
u[:,-1] = Vinf # free stream at top

bl_thckns = np.zeros(Nx)

for i in range(0,Nx-1):
    for j in range(0,Ny-2):
        dudy = 0.5 * (u[i, j + 2] - u[i, j])
        d2udy2 = (u[i, j + 2] - 2 * u[i, j + 1] + u[i, j]) / u[i, j + 1]
        u[i+1,j+1] = u[i,j+1] + nu*d2udy2*(dx/dy**2) - dudy*(v[i,j+1]/u[i,j+1])*dx/dy
        v[i+1,j+1] = v[i+1,j] - 0.5*(u[i+1,j+1]-u[i,j+1]+u[i+1,j]-u[i,j])*dy/dx


for t in range(0,Nx):
    for k in range(0,Ny):
        if u[t,k] > 0.99*Vinf:
            bl_thckns[t] = k*dy
            break

x = np.linspace(0, L, Nx)
y = np.linspace(0, H, Ny)


d1 = np.zeros(Nx)
d2 = np.zeros(Nx)
dudytau = np.zeros(Nx)
tw = np.zeros(Nx)
cf = np.zeros(Nx)
for w in range(Nx):
    d1[w] = np.sum((1-u[w,:]/Vinf)*dy)
    d2[w] = np.sum((u[w,:] / Vinf) * (1 - u[w, :] / Vinf) * dy)
    dudytau[w] = (u[w,1]-u[w,0])/dy
    tw[w] = nu*rho*dudytau[w]
    cf[w] = tw[w]/(0.5*rho*Vinf**2)

plt.plot(x, bl_thckns)
plt.xlabel("x")
plt.ylabel("Boundary Layer Thickness")
plt.show()


def run_explicit():
    return x, y, bl_thckns, d1, d2, cf, u, v