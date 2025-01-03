import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from scipy.optimize import fsolve


# def conservation_check(rho, V, T, A, t, dx):
#     mass = np.sum(rho[t, :] * A) * dx
#     momentum = np.sum(rho[t, :] * V[t, :] * A) * dx
#     energy = np.sum(rho[t, :] * T[t, :] * A) * dx
#     return mass, momentum, energy


gamma = 1.4
L = 3
Nx = 31
Nt = 5000
x = np.linspace(0,L,Nx) # x/L
dx = L/(Nx-1)
C = 0.5

p0 = 1
pe = 0.93


# initials
rho_i = 1 - 0.023*x # rho/rho_0
T_i = 1 - 0.009333*x # T/T_0
V_i = 0.05+0.11*x # V/a_0

# nozzle geometry
limit1 = np.where(x == 1.5)[0][0]
A = np.zeros(len(x))
for i in range(0,limit1):
    A[i] = 1+2.2*(x[i] - 1.5)**2 # A/A*
for i in range(limit1, len(x)):
    A[i] = 1 + 0.2223 * (x[i] - 1.5) ** 2

# A = A_/min(A_)

dt = np.min(C*dx/(T_i**0.5+V_i))
tss = (Nt-1)*dt
time = np.linspace(0, tss, int(Nt))
# print(len(time))


rho = np.zeros((len(time),len(x)))
V = np.zeros((len(time),len(x)))
T = np.zeros((len(time),len(x)))
p = np.zeros((len(time),len(x)))
M = np.zeros((len(time),len(x)))
drhodt = np.zeros((len(time),len(x)))
dVdt = np.zeros((len(time),len(x)))
dTdt = np.zeros((len(time),len(x)))


rho_est = np.zeros((len(time),len(x)))
V_est = np.zeros((len(time),len(x)))
T_est = np.zeros((len(time),len(x)))
drhodt_est = np.zeros((len(time),len(x)))
dVdt_est = np.zeros((len(time),len(x)))
dTdt_est = np.zeros((len(time),len(x)))
drhodt_av = np.zeros((len(time),len(x)))
dVdt_av = np.zeros((len(time),len(x)))
dTdt_av = np.zeros((len(time),len(x)))

rho_history = []
p_history = []
T_history = []
M_history = []
drhodt_av_history = []
dVdt_av_history = []

mass_flow = np.zeros((len(time),len(x)))

mid = (Nx-1)/2
rho[0,:] = rho_i
V[0,:] = V_i
T[0,:] = T_i
mass_flow[0,:] = rho[0,:]*V[0,:]*A[:]

# # Analytical calculations
# def Mach_eq(M_an, A_):
#     return (A_) ** 2 - (1 / M_an ** 2) * (
#                 ((2 / (gamma + 1)) * (1 + ((gamma - 1) / 2) * M_an ** 2)) ** ((gamma + 1) / (gamma - 1)))
#
#
# Mtot = np.zeros(Nx)
#
# for i in range(0, Nx):
#     if i < (Nx) / 2:
#         init_guess = 0.2
#     else:
#         init_guess = 2
#     M_an = fsolve(Mach_eq, init_guess, args=A_[i])
#
#     if M_an < 0:
#         M_an = -M_an
#     Mtot[i] = M_an
#
# p_an = (1 + (gamma - 1) / 2 * Mtot ** 2) ** (-gamma / (gamma - 1))
# rho_an = (1 + (gamma - 1) / 2 * Mtot ** 2) ** (-1 / (gamma - 1))
# T_an = (1 + (gamma - 1) / 2 * Mtot ** 2) ** -1
#
#

p[:,0] = p0
p[:,-1] = pe


for t in range(0,Nt-1):

    for i in range(0,Nx-1):
        # derivatives
        drhodt[t, i] = - rho[t,i] * (V[t,i+1]-V[t,i])/dx - rho[t,i]*V[t,i]*(np.log(A[i+1])-np.log(A[i]))/dx - V[t,i]*(rho[t,i+1]-rho[t,i])/dx
        dVdt[t,i] = -V[t,i]*(V[t,i+1]-V[t,i])/dx - 1/gamma *((T[t,i+1]-T[t,i])/dx + T[t,i]/rho[t,i]*(rho[t,i+1]-rho[t,i])/dx)
        dTdt[t, i] = -V[t,i]*(T[t,i+1]-T[t,i])/dx - (gamma-1)*T[t,i]*((V[t,i+1]-V[t,i])/dx + V[t,i]*(np.log(A[i+1]) - np.log(A[i]))/dx)

    for i in range(0,Nx):
        # estimators
        rho_est[t+1,i] = rho[t, i] + drhodt[t,i]*dt
        V_est[t + 1, i] = V[t, i] + dVdt[t, i] * dt
        T_est[t + 1, i] = T[t, i] + dTdt[t, i] * dt

    for i in range(1,Nx):
        # der estimators
        drhodt_est[t+1, i] = - rho_est[t+1, i] * (V_est[t+1, i] - V_est[t+1, i-1]) / dx - rho_est[t+1, i] * V_est[t+1, i] * (np.log(A[i]) - np.log(A[i-1])) / dx - V_est[t+1, i] * (rho_est[t+1, i] - rho_est[t+1, i-1]) / dx
        dVdt_est[t+1, i] = - V_est[t + 1, i] * (V_est[t + 1, i] - V_est[t + 1, i-1]) / dx - 1/gamma* ((T_est[t + 1, i]-T_est[t+1,i-1])/dx + T_est[t + 1, i]/rho_est[t+1,i] * (rho_est[t + 1, i] - rho_est[t + 1, i-1]) / dx)
        dTdt_est[t+1, i] = - V_est[t + 1, i] * (T_est[t + 1, i] - T_est[t + 1, i-1]) / dx - (gamma-1) * T_est[t + 1, i] * ((V_est[t + 1, i] - V_est[t + 1, i-1]) / dx + V_est[t+1,i]*(np.log(A[i]) - np.log(A[i-1])) / dx)

    for i in range(1,Nx-1):
        drhodt_av[t, i] = (drhodt_est[t + 1, i] + drhodt[t, i]) / 2
        dVdt_av[t, i] = (dVdt_est[t + 1, i] + dVdt[t, i])*0.5
        dTdt_av[t, i] = (dTdt_est[t + 1, i] + dTdt[t, i]) / 2
    drhodt_av[t,0] = drhodt_est[t+1,0]
    dVdt_av[t,0] = dVdt_est[t+1,0]
    dTdt_av[t,0] = dTdt_est[t+1,0]
    drhodt_av[t,-1] = drhodt_est[t,-1]
    dVdt_av[t,-1] = dVdt_est[t,-1]
    dTdt_av[t,-1] = dTdt_est[t,-1]

    for i in range(0,Nx):
        # corrector
        rho[t+1,i] = rho[t,i] + drhodt_av[t,i]*dt
        V[t + 1, i] = V[t, i] + dVdt_av[t,i] * dt
        T[t + 1, i] = T[t, i] + dTdt_av[t,i] * dt

    # edge values
    rho[t + 1, 0] = 1
    T[t + 1, 0] = 1
    V[t + 1, 0] = 2 * V[t + 1, 1] - V[t + 1, 2]
    V[t + 1, -1] = 2 * V[t + 1, -2] - V[t + 1, -3]
    T[t + 1, -1] = 2 * T[t + 1, -2] - T[t + 1, -3]
    rho[t + 1, -1] = pe / T[t + 1, -1]
    # pressure and Mach
    p[t + 1,:] = rho[t + 1,:] * T[t + 1,:]
    M[t+1,:] = V[t+1,:] / (T[t + 1,:]**0.5)



    # if t % 100 == 0:  # Check conservation every 100 time steps
    #     mass, momentum, energy = conservation_check(rho, V, T, A, t, dx)
    #     print(f"Time step {t}, Mass = {mass:.5f}, Momentum = {momentum:.5f}, Energy = {energy:.5f}")



    # plot vectors
    rho_history.append(rho[t + 1, int(mid)])
    p_history.append(p[t+1,int(mid)])
    T_history.append(T[t+1,int(mid)])
    M_history.append(M[t+1,int(mid)])

    drhodt_av_history.append(np.abs(drhodt_av[t,int(mid)]))
    dVdt_av_history.append(np.abs(dVdt_av[t,int(mid)]))

    # for i in range(0,Nx):
    #     mass_flow[t+1,i] = rho[t+1,i]*V[t+1,i]*A[i]

print(rho[4999,15], V[4999,15], T[4999,15])

# ########## RESULTS ###########
#
# # Tab. 7.3
# print(np.round(np.array((x.T, A.T, rho[1399,:].T, V[1399,:].T, T[1399,:].T, p[1399,:].T, M[1399,:].T, mass_flow[1399,:].T)),3))
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







