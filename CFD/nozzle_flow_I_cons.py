import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve


C = 0.9  # Courant
gamma = 1.4
L = 3  # m
Nx = 41
dx = L / (Nx-1)
Nt = 1900
x = np.linspace(0, L, Nx)
A_ = 1 + 2.2 * (x - 1.5) ** 2
A = A_/min(A_)

# initials
limit1 = np.where(x >= 0.5)[0][0]
limit2 = np.where(x >= 1.5)[0][0]

rho = np.zeros(len(x))
T = np.zeros(len(x))

for i in range(0,limit1):
    rho[i] = 1
    T[i] = 1
for i in range(limit1,limit2):
    rho[i] = 1-0.366*(x[i]-0.5)
    T[i] = 1 - 0.167 * (x[i] - 0.5)
for i in range(limit2,Nx):
    rho[i] = 0.634-0.3879*(x[i]-1.5)
    T[i] = 0.833-0.3507*(x[i]-1.5)


V = 0.59 / (rho * A)
p = rho*T
U1 = rho * A
U2 = rho * A * V
U3 = rho * (T / (gamma - 1) + gamma / 2 * V ** 2) * A

F1 = U2
F2 = U2 ** 2 / U1 + (gamma - 1) / gamma * (U3 - gamma / 2 * U2 ** 2 / U1)
F3 = gamma * U2 * U3 / U1 - gamma * (gamma - 1) / 2 * U2 ** 3 / U1 ** 2

J2 = np.zeros(len(x) - 1)
dU1dt = np.zeros(len(x) - 1)
dU2dt = np.zeros(len(x) - 1)
dU3dt = np.zeros(len(x) - 1)
dt = min(C * dx / (T[:] ** 0.5 + V[:]))
rho_est = np.zeros(len(x) - 1)
T_est = np.zeros(len(x) - 1)

U1_est = np.zeros(len(x) - 1)
U2_est = np.zeros(len(x) - 1)
U3_est = np.zeros(len(x) - 1)
F1_est = np.zeros(len(x) - 1)
F2_est = np.zeros(len(x) - 1)
F3_est = np.zeros(len(x) - 1)

mass = {}
mass[0] = rho*A*V
pressure = {}
mid = int((Nx-1)/2)


# Analytical calculations
def M_solve(M_, A_):
    return (A_) ** 2 - (1 / M_ ** 2) * (
                ((2 / (gamma + 1)) * (1 + ((gamma - 1) / 2) * M_ ** 2)) ** ((gamma + 1) / (gamma - 1)))


M_an = np.zeros(Nx)

for i in range(Nx):
    guess = 0.2 if x[i] < 1.5 else 2.0
    M_ = fsolve(M_solve, guess, args=(A_[i],))[0]
    M_an[i] = abs(M_)

p_an = (1 + (gamma - 1) / 2 * M_an ** 2) ** (-gamma / (gamma - 1))
rho_an = (1 + (gamma - 1) / 2 * M_an ** 2) ** (-1 / (gamma - 1))
T_an = (1 + (gamma - 1) / 2 * M_an ** 2) ** -1




rho_history = []
T_history = []
p_history = []
M_history = []
for j in range(Nt):
    for i in range(len(x) - 1):
        # Predictor Step
        J2[i] = (rho[i] * T[i]) / gamma * ((A[i+1] - A[i]) / dx)
        dU1dt[i] = -(F1[i + 1] - F1[i]) / dx
        dU2dt[i] = -(F2[i + 1] - F2[i]) / dx + J2[i]
        dU3dt[i] = -(F3[i + 1] - F3[i]) / dx
        U1_est[i] = U1[i] + dU1dt[i] * dt
        U2_est[i] = U2[i] + dU2dt[i] * dt
        U3_est[i] = U3[i] + dU3dt[i] * dt
        rho_est[i] = U1_est[i] / A[i]
        T_est[i] = (gamma - 1) * (U3_est[i] / U1_est[i] - gamma / 2 * (U2_est[i] / U1_est[i]) ** 2)
        F1_est[i] = U2_est[i]
        F2_est[i] = U2_est[i] ** 2 / U1_est[i] + (gamma - 1) / gamma * (U3_est[i] - gamma / 2 * U2_est[i] ** 2 / U1_est[i])
        F3_est[i] = gamma * U2_est[i] * U3_est[i] / U1_est[i] - gamma * (gamma - 1) / 2 * (U2_est[i] / U1_est[i]) ** 2 * U2_est[i]

        # Corrector Step
    for i in range(1, len(x) - 1):
        dU1dt_est = -(F1_est[i] - F1_est[i - 1]) / dx
        dU2dt_est = -(F2_est[i] - F2_est[i - 1]) / dx + 1 / gamma * rho_est[i] * T_est[i] * (A[i] - A[i - 1]) / dx
        dU3dt_est = -(F3_est[i] - F3_est[i - 1]) / dx
        dU1_av = 0.5 * (dU1dt[i] + dU1dt_est)
        dU2_av = 0.5 * (dU2dt[i] + dU2dt_est)
        dU3_av = 0.5 * (dU3dt[i] + dU3dt_est)

        U1[i] = U1[i] + dU1_av * dt
        U2[i] = U2[i] + dU2_av * dt
        U3[i] = U3[i] + dU3_av * dt

        rho[i] = U1[i] / A[i]
        V[i] = U2[i] / U1[i]
        T[i] = (gamma - 1) * (U3[i] / U1[i] - gamma / 2 * V[i] ** 2)

        # Boundary Conditions
        U1[0] = A[0]
        U2[0] = 2 * U2[1] - U2[2]
        V[0] = U2[0] / U1[0]
        U3[0] = U1[0] * (T[0] / (gamma - 1) + gamma / 2 * V[0] ** 2)
        rho[0] = 1
        T[0] = 1


        U1[-1] = 2 * U1[-2] - U1[-3]
        U2[-1] = 2 * U2[-2] - U2[-3]
        U3[-1] = 2 * U3[-2] - U3[-3]
        rho[-1] = U1[-1] / A[-1]
        V[-1] = U2[-1]/U1[-1]
        T[-1] = (gamma - 1) * (U3[-1] / U1[-1] - gamma / 2 * V[-1] ** 2)

    mass_flow = rho * V * A

    F1 = U2
    F2 = U2 ** 2 / U1 + (gamma - 1) / gamma * (U3 - gamma / 2 * U2 ** 2 / U1)
    F3 = gamma * U2 * U3 / U1 - gamma * (gamma - 1) / 2 * U2 ** 3 / U1 ** 2

    p = rho*T
    M = V/(T**0.5)

    if j == 0 or j == 49 or j == 99 or j == 149 or j == 199 or j == 699:
        mass[j] = mass_flow

    rho_history.append(rho[mid])
    T_history.append(T[mid])
    p_history.append(p[mid])
    M_history.append(M[mid])



# Tab. 7.3
print(x.T[:10], A.T[:10], rho[:10].T, V[:10].T, T[:10].T, p[:10].T, M[:10].T, mass_flow[:10].T)

# Tab. 7.4
print(np.round(np.array((x.T[:10], A.T[:10], rho[:10].T, rho_an[:10], np.abs(rho[:10]-rho_an[:10])/rho[:10]*100, M[:10].T, M_an[:10].T, np.abs(M[:10]-M_an[:10])/M[:10]*100)),3))

# Tab. 7.5
# results can be found on Tab. 7.6 if the grid points are changed

# # # Tab. 7.6
print(f"Density numerical = {rho[int(mid)]}, Density analytical = {rho_an[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")
print(f"Temperature numerical = {T[int(mid)]}, Temperature analytical = {T_an[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")
print(f"Pressure numerical = {p[int(mid)]}, Pressure analytical = {p_an[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")
print(f"Mach numerical = {M[int(mid)]}, Mach analytical = {M_an[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")



#
# # Fig. 7.9
plt.figure(figsize=(10, 6))
plt.plot(range(1, len(rho_history) + 1), rho_history, label="At Nozzle Throat")
plt.xlabel("Number of Time Steps")
plt.ylabel(r"$\rho / \rho_0$")
plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(range(1, len(T_history) + 1), T_history, label="At Nozzle Throat")
plt.xlabel("Number of Time Steps")
plt.ylabel(r"$T / T_0$")
plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(range(1, len(p_history) + 1), p_history, label="At Nozzle Throat")
plt.xlabel("Number of Time Steps")
plt.ylabel(r"$p / p_0$")
plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(range(1, len(M_history) + 1), M_history, label="At Nozzle Throat")
plt.xlabel("Number of Time Steps")
plt.ylabel(r"$M$")
plt.legend()
plt.grid()
plt.show()
#
# Fig. 7.10
# plt.figure(figsize=(10, 6))
# plt.plot(range(1, len(drdt_av_history) + 1), drdt_av_history, label=r"$|(\frac{d\rho}{dt})_{avg}|$")
# plt.plot(range(1, len(dVdt_av_history) + 1), dVdt_av_history, label=r"$|(\frac{dV}{dt})_{avg}|$")
# plt.xlabel("Number of Time Steps")
# plt.ylabel("Residual")
# plt.title(f"Dimensionless time derivatives at Grid Point {int(mid+1)} Over Time")
# plt.legend()
# plt.grid()
# plt.show()



# # Fig. 7.11
plt.figure(figsize=(10, 6))
plt.plot(x,mass[0], label=r"$0\Delta t$")
plt.plot(x,mass[49], label=r"$50\Delta t$")
plt.plot(x,mass[99], label=r"$100\Delta t$")
plt.plot(x,mass[149], label=r"$150\Delta t$")
plt.plot(x,mass[199], label=r"$150\Delta t$")
plt.plot(x,mass[699], label=r"$700\Delta t$")
plt.xlabel("Nondimensionless distance through nozzle (x)")
plt.ylabel("Nondimensionless mass flow")
plt.legend()
plt.grid()
plt.show()
#
# # Fig. 7.12
fig, ax1 = plt.subplots()
ax1.plot(x, rho, label="Numerical results")
ax1.plot(x[::3], rho_an[::3], 'o', label="Analytical results")
ax1.set_xlabel("Nondimensionless distance through nozzle (x)")
ax1.set_ylabel("Nondimensionless density")
ax2 = ax1.twinx()
ax2.plot(x, M, 'g', label="Numerical results")
ax2.plot(x[::3], M_an[::3],'o', label="Analytical results")
ax2.set_ylabel("Mach number")
plt.legend()
plt.grid()
plt.show()
