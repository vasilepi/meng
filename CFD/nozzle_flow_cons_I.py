import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from scipy.optimize import fsolve


def conservation_check(rho, V, T, A, t, dx):
    mass = np.sum(rho[t, :] * A) * dx
    momentum = np.sum(rho[t, :] * V[t, :] * A) * dx
    energy = np.sum(rho[t, :] * T[t, :] * A) * dx
    return mass, momentum, energy


gamma = 1.4
L = 3
Nx = 31
Nt = 202
x = np.linspace(0,L,Nx) # x/L
dx = L/(Nx-1)
C = 0.5

# initials
limit1 = np.where(x == 0.5)[0][0]
limit2 = np.where(x == 1.5)[0][0]

rho_i = np.zeros(len(x))
T_i = np.zeros(len(x))
V_i = np.zeros(len(x))
U1_i = np.zeros(len(x))
U2_i = np.ones(len(x))*0.59
U3_i = np.zeros(len(x))

for i in range(0,limit1):
    rho_i[i] = 1
    T_i[i] = 1
for i in range(limit1,limit2):
    rho_i[i] = 1-0.366*(x[i]-0.5)
    T_i[i] = 1 - 0.167 * (x[i] - 0.5)
for i in range(limit2,Nx):
    rho_i[i] = 0.634-0.3879*(x[i]-1.5)
    T_i[i] = 0.833-0.3507*(x[i]-1.5)


# nozzle geometry
A_ = 1+2.2*(x - 1.5)**2 # A/A*
A = A_/min(A_)


U1_i = rho_i*A
V_i = U2_i/U1_i
U3_i = rho_i*(T_i/(gamma-1) + 0.5*gamma*V_i**2)*A



# dt = np.min(C*dx/(T_i**0.5+V_i))
dt = 0.0267
tss = (Nt-1)*dt
time = np.linspace(0, tss, int(Nt))



rho = np.zeros((len(time),len(x)))
V = np.zeros((len(time),len(x)))
T = np.zeros((len(time),len(x)))
p = np.zeros((len(time),len(x)))
M = np.zeros((len(time),len(x)))
U1 = np.zeros((len(time),len(x)))
U2 = np.zeros((len(time),len(x)))
U3 = np.zeros((len(time),len(x)))
F1 = np.zeros((len(time),len(x)))
F2 = np.zeros((len(time),len(x)))
F3 = np.zeros((len(time),len(x)))
J2 = np.zeros((len(time),len(x)))
dU1dt = np.zeros((len(time),len(x)))
dU2dt = np.zeros((len(time),len(x)))
dU3dt = np.zeros((len(time),len(x)))
#
#
rho_est = np.zeros((len(time),len(x)))
U1_est = np.zeros((len(time),len(x)))
U2_est = np.zeros((len(time),len(x)))
U3_est = np.zeros((len(time),len(x)))
F1_est = np.zeros((len(time),len(x)))
F2_est = np.zeros((len(time),len(x)))
F3_est = np.zeros((len(time),len(x)))
T_est = np.zeros((len(time),len(x)))
J2_est = np.zeros((len(time),len(x)))
dU1dt_est = np.zeros((len(time),len(x)))
dU2dt_est = np.zeros((len(time),len(x)))
dU3dt_est = np.zeros((len(time),len(x)))
dU1dt_av = np.zeros((len(time),len(x)))
dU2dt_av = np.zeros((len(time),len(x)))
dU3dt_av = np.zeros((len(time),len(x)))


rho[0,:] = rho_i
V[0,:] = V_i
T[0,:] = T_i
U1[0,:] = U1_i
U2[0,:] = U2_i
U3[0,:] = U3_i

F1[0,:] = U2[0,:]
F2[0,:] = (U2[0,:]**2/U1[0,:] + (gamma-1)/gamma * (U3[0,:]-0.5*gamma*U2[0,:]**2/U1[0,:]))
F3[0,:] = U2[0,:]*U3[0,:]*gamma/U1[0,:] - gamma*(gamma-1)*0.5*U2[0,:]**3/(U1[0,:]**2)
for i in range(0,Nx-1):
    J2[0,i] = 1/gamma * rho[0,i]*T[0,i]*(A[i+1]-A[i])/dx



# rho[:,0] = 1
# T[:,0] = 1
# U1[:,0] = 1
# U2[0,0] = 2*U2[0,1]-U2[0,2]
# U3[0,0] = U1[0,0]*(T[0,0]/(gamma-1)+0.5*gamma*V[0,0]*2)


for t in range(Nt-1):
# for t in range(0,1):

    for i in range(0,Nx-1):
        # derivatives
        dU1dt[t, i] = - (F1[t,i+1]-F1[t,i])/dx
        dU2dt[t,i] = - (F2[t,i+1]-F2[t,i])/dx + J2[t,i]
        dU3dt[t, i] = - (F3[t,i+1]-F3[t,i])/dx

    for i in range(0,Nx-1):
        # estimators
        U1_est[t+1,i] = U1[t, i] + dU1dt[t,i]*dt
        U2_est[t + 1, i] = U2[t, i] + dU2dt[t, i] * dt
        U3_est[t + 1, i] = U3[t, i] + dU3dt[t, i] * dt


        rho_est[t+1,i] = U1_est[t+1,i]/A[i]
        T_est[t + 1, i] = (gamma-1)*(U3_est[t+1,i]/U1_est[t + 1, i]-0.5*gamma*(U2_est[t + 1, i]/U1_est[t + 1, i])**2)


        F1_est[t+1, i] = U2_est[t+1, i]
        F2_est[t+1,i] = (U2_est[t+1, i]**2 / U1_est[t+1, i] + (gamma - 1) / gamma * (U3_est[t+1,i] - 0.5 * gamma * U2_est[t+1,i]**2 / U1_est[t+1, i]))
        F3_est[t+1, i] = U2_est[t+1,i] * U3_est[t+1,i] * gamma / U1_est[t+1,i] - gamma * (gamma - 1) * 0.5 * U2_est[t+1,i]**3 / (U1_est[t+1,i]**2)

    for i in range(1,Nx):
        # J2_est[t+1,i] = 1 / gamma * rho_est[t + 1, i] * T_est[t + 1, i] * (A[i+1] - A[i]) / dx
        J2_est[t + 1, i] = (gamma-1) / gamma * (U3_est[t + 1, i] * 0.5*gamma*U2_est[t + 1, i]**2/U1_est[t+1,i]) * (np.log(A[i]) - np.log(A[i-1])) / dx
    for i in range(1,Nx):
        # der estimators
        dU1dt_est[t+1, i] = - (F1_est[t+1,i]-F1_est[t+1,i-1])/dx
        dU2dt_est[t+1, i] = - (F2[t+1,i]-F2[t+1,i-1])/dx + J2_est[t+1,i]
        dU3dt_est[t+1, i] = - (F3_est[t+1,i]-F3_est[t+1,i-1])/dx

    for i in range(1,Nx-1):
        dU1dt_av[t, i] = (dU1dt_est[t + 1, i] + dU1dt[t, i]) / 2
        dU2dt_av[t, i] = (dU2dt_est[t + 1, i] + dU2dt[t, i]) / 2
        dU3dt_av[t, i] = (dU3dt_est[t + 1, i] + dU3dt[t, i]) / 2

    # dU1dt_av[t,0] = dU1dt_est[t+1,0]
    # dU2dt_av[t,0] = dU2dt_est[t+1,0]
    # dU3dt_av[t,0] = dU3dt_est[t+1,0]
    # dU1dt_av[t,-1] = dU1dt_est[t+1,-1]
    # dU2dt_av[t,-1] = dU2dt_est[t+1,-1]
    # dU3dt_av[t,-1] = dU3dt_est[t+1,-1]

    for i in range(1,Nx-1):
        # corrector
        U1[t+1,i] = U1[t,i] + dU1dt_av[t,i]*dt
        U2[t + 1, i] = U2[t, i] + dU2dt_av[t,i]*dt
        U3[t + 1, i] = U3[t, i] + dU3dt_av[t,i]*dt



    for i in range(1,Nx-1):
        rho[t+1,i] = U1[t+1,i]/A[i]
    for i in range(1,Nx-1):
        V[t+1,i] = U2[t+1,i]/U1[t+1,i]
    for i in range(1,Nx-1):
        T[t+1,i] = (gamma-1)*(U3[t+1,i]/U1[t + 1, i]-0.5*gamma*V[t+1,i]**2)


    rho[t+1,0] = 1
    T[t+1,0] = 1
    U1[t+1,0] = A[0]
    U2[t+1,0] = 2*U2[t+1,1]-U2[t+1,2]
    U3[t+1,0] = U1[t+1,0]*(T[t+1,0]/(gamma-1)+0.5*gamma*V[t+1,0]*2)
    F1[t+1,0] = U2[t+1,0]
    F2[t+1,0] = (U2[t+1,0]**2/U1[t+1,0] + (gamma-1)/gamma * (U3[t+1,0]-0.5*gamma*U2[t+1,0]**2/U1[t+1,0]))
    F3[t+1,0] = U2[t+1,0]*U3[t+1,0]*gamma/U1[t+1,0] - gamma*(gamma-1)*0.5*U2[t+1,0]**3/(U1[t+1,0]**2)


    U1[t+1,-1] = 2*U1[t+1,-2]-U1[t+1,-3]
    U2[t+1,-1] = 2*U2[t+1,-2]-U2[t+1,-3]
    U3[t+1,-1] = 2*U3[t+1,-2]-U3[t+1,-3]
    F1[t+1,-1] = U2[t+1,-1]
    F2[t+1,-1] = (U2[t+1,-1]**2/U1[t+1,-1] + (gamma-1)/gamma * (U3[t+1,-1]-0.5*gamma*U2[t+1,-1]**2/U1[t+1,-1]))
    F3[t+1,-1] = U2[t+1,-1]*U3[t+1,-1]*gamma/U1[t+1,-1] - gamma*(gamma-1)*0.5*U2[t+1,-1]**3/(U1[t+1,-1]**2)

    for i in range(0,Nx-1):
        J2[t+1,i] = (gamma-1) / gamma * (U3[t + 1, i] * 0.5*gamma*U2[t + 1, i]**2/U1[t+1,i]) * (np.log(A[i + 1]) - np.log(A[i])) / dx

    if t % 100 == 0:  # Check conservation every 100 time steps
        mass, momentum, energy = conservation_check(rho, V, T, A, t, dx)
        print(f"Time step {t}, Mass = {mass:.5f}, Momentum = {momentum:.5f}, Energy = {energy:.5f}")







print(rho[1,15], T[1,15], V[1,15], U1[1,15], U2[1,15], U3[1,15], U1_est[1,15], F1_est[1,14], U2_est[1,15], F2_est[1,14], U3_est[1,15], F3_est[1,14])
# ########## RESULTS ###########
#
# # Tab. 7.3
# print(np.round(np.array((x.T, A.T, rho[1,:].T, V[1,:].T, T[1,:].T)),3))
#
# # Tab. 7.4
# print(np.round(np.array((x.T, A.T, rho[1399,:].T, rho_an[:], np.abs(rho[1399,:]-rho_an[:])/rho[1399,:]*100, M[1399,:].T, Mtot[:].T, np.abs(M[1399,:]-Mtot[:])/M[1399,:]*100)),3))
#
# # Tab. 7.5
# # results can be found on Tab. 7.6 if the grid points are changed
#
# # Tab. 7.6
# print(f"Density numerical = {rho[1399,int(mid)]}, Density analytical = {rho_an[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")
# print(f"Temperature numerical = {T[1399,int(mid)]}, Temperature analytical = {T_an[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")
# print(f"Pressure numerical = {p[1399,int(mid)]}, Pressure analytical = {p_an[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")
# print(f"Mach numerical = {M[1399,int(mid)]}, Mach analytical = {Mtot[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")
#
#
#
#
# # Fig. 7.9
# plt.figure(figsize=(10, 6))
# plt.plot(range(1, len(rho_history) + 1), rho_history, label="Numerical Solution")
# plt.xlabel("Number of Time Steps")
# plt.ylabel(r"$\rho / \rho_0$")
# plt.title(f"Density at Grid Point {int(mid+1)} Over Time")
# plt.legend()
# plt.grid()
# plt.show()
#
# plt.figure(figsize=(10, 6))
# plt.plot(range(1, len(T_history) + 1), T_history, label="Numerical Solution")
# plt.xlabel("Number of Time Steps")
# plt.ylabel(r"$T / T_0$")
# plt.title(f"Temperature at Grid Point {int(mid+1)} Over Time")
# plt.legend()
# plt.grid()
# plt.show()
#
# plt.figure(figsize=(10, 6))
# plt.plot(range(1, len(p_history) + 1), p_history, label="Numerical Solution")
# plt.xlabel("Number of Time Steps")
# plt.ylabel(r"$p / p_0$")
# plt.title(f"Pressure at Grid Point {int(mid+1)} Over Time")
# plt.legend()
# plt.grid()
# plt.show()
#
# plt.figure(figsize=(10, 6))
# plt.plot(range(1, len(M_history) + 1), M_history, label="Numerical Solution")
# plt.xlabel("Number of Time Steps")
# plt.ylabel(r"$M$")
# plt.title(f"Mach number at Grid Point {int(mid+1)} Over Time")
# plt.legend()
# plt.grid()
# plt.show()
#
# # Fig. 7.10
# plt.figure(figsize=(10, 6))
# plt.plot(range(1, len(drhodt_av_history) + 1), drhodt_av_history, label=r"$|(\frac{d\rho}{dt})_{avg}|$")
# plt.plot(range(1, len(dVdt_av_history) + 1), dVdt_av_history, label=r"$|(\frac{dV}{dt})_{avg}|$")
# plt.xlabel("Number of Time Steps")
# plt.ylabel("Residual")
# plt.title(f"Dimensionless time derivatives at Grid Point {int(mid+1)} Over Time")
# plt.legend()
# plt.grid()
# plt.show()
#
# # Fig. 7.11
# plt.figure(figsize=(10, 6))
# plt.plot(x,mass_flow[0,:], label=r"$0\Delta t$")
# plt.plot(x,mass_flow[50,:], label=r"$50\Delta t$")
# plt.plot(x,mass_flow[100,:], label=r"$100\Delta t$")
# plt.plot(x,mass_flow[150,:], label=r"$150\Delta t$")
# plt.plot(x,mass_flow[200,:], label=r"$200\Delta t$")
# plt.plot(x,mass_flow[700,:], label=r"$700\Delta t$")
# plt.xlabel("Nondimensionless distance through nozzle (x)")
# plt.ylabel("Nondimensionless mass flow")
# plt.title("Mass flow distributions")
# plt.legend()
# plt.grid()
# plt.show()
#
# # Fig. 7.12
# plt.figure(figsize=(10,6))
# fig, ax1 = plt.subplots()
# ax1.plot(x, rho[1399,:], label="Numerical results")
# ax1.plot(x[::3], rho_an[::3], 'o', label="Analytical results")
# ax1.set_xlabel("Nondimensionless distance through nozzle (x)")
# ax1.set_ylabel("Nondimensionless density")
# ax2 = ax1.twinx()
# ax2.plot(x, M[1399,:], 'g', label="Numerical results")
# ax2.plot(x[::3], Mtot[::3],'o', label="Analytical results")
# ax2.set_ylabel("Mach number")
# plt.title("Steady-state distributions over nozzle distance")
# plt.legend()
# plt.grid()
# plt.show()







