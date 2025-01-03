import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import fsolve


C = 0.5  # Courant
g = 1.4
L = 3  # m
Nx = 31
dx = L / (Nx-1)
Nt = 1400
x = np.linspace(0, L, Nx)
A_ = 1 + 2.2 * (x - 1.5) ** 2
A = A_/min(A_)
r = np.ones(len(x))
T = np.ones(len(x))
r[int(0.5 / dx):int(1.5 / dx) + 1] = 1 - 0.366 * (x[int(0.5 / dx):int(1.5 / dx) + 1] - 0.5)
T[int(0.5 / dx):int(1.5 / dx) + 1] = 1 - 0.167 * (x[int(0.5 / dx):int(1.5 / dx) + 1] - 0.5)
r[int(1.5 / dx):int(L / dx) + 1] = 0.634 - 0.3879 * (x[int(1.5 / dx):int((L + dx) / dx)] - 1.5)
T[int(1.5 / dx):int(L / dx) + 1] = 0.833 - 0.3507 * (x[int(1.5 / dx):int((L + dx) / dx)] - 1.5)
V = 0.59 / (r * A)
U1 = r * A
U2 = r * A * V
U3 = r * (T / (g - 1) + g / 2 * V ** 2) * A

F1 = U2
F2 = U2 ** 2 / U1 + (g - 1) / g * (U3 - g / 2 * U2 ** 2 / U1)
F3 = g * U2 * U3 / U1 - g * (g - 1) / 2 * U2 ** 3 / U1 ** 2

J2 = np.zeros(len(x) - 1)
dU1 = np.zeros(len(x) - 1)
dU2 = np.zeros(len(x) - 1)
dU3 = np.zeros(len(x) - 1)
Dt = min(C * dx / (T[:] ** 0.5 + V[:]))
r_ = np.zeros(len(x) - 1)
T_ = np.zeros(len(x) - 1)

U1_ = np.zeros(len(x) - 1)
U2_ = np.zeros(len(x) - 1)
U3_ = np.zeros(len(x) - 1)
F1_ = np.zeros(len(x) - 1)
F2_ = np.zeros(len(x) - 1)
F3_ = np.zeros(len(x) - 1)

mass = {}
mass[0] = r*A*V
pressure = {}
mid = (Nx-1)/2
# Analytical calculations
def Mach_eq(M_an, A_):
    return (A_) ** 2 - (1 / M_an ** 2) * (
                ((2 / (g + 1)) * (1 + ((g - 1) / 2) * M_an ** 2)) ** ((g + 1) / (g - 1)))


Mtot = np.zeros(Nx)

for i in range(0, Nx):
    if i < (Nx) / 2:
        init_guess = 0.2
    else:
        init_guess = 2
    M_an = fsolve(Mach_eq, init_guess, args=A_[i])

    if M_an < 0:
        M_an = -M_an
    Mtot[i] = M_an

p_an = (1 + (g - 1) / 2 * Mtot ** 2) ** (-g / (g - 1))
r_an = (1 + (g - 1) / 2 * Mtot ** 2) ** (-1 / (g - 1))
T_an = (1 + (g - 1) / 2 * Mtot ** 2) ** -1

r_history = []
T_history = []
p_history = []
M_history = []
for j in range(Nt):
    for i in range(len(x) - 1):
        # Predictor Step
        # J2[i] = (g-1)/g * (U3[i] - g/2 * U2[i]**2/U1[i]) * (A[i]-A[i-1])/dx
        J2[i] = 1 / g * r[i] * T[i] * (A[i + 1] - A[i]) / dx
        dU1[i] = -(F1[i + 1] - F1[i]) / dx
        dU2[i] = -(F2[i + 1] - F2[i]) / dx + J2[i]
        dU3[i] = -(F3[i + 1] - F3[i]) / dx
        U1_[i] = U1[i] + dU1[i] * Dt
        U2_[i] = U2[i] + dU2[i] * Dt
        U3_[i] = U3[i] + dU3[i] * Dt
        r_[i] = U1_[i] / A[i]
        T_[i] = (g - 1) * (U3_[i] / U1_[i] - g / 2 * (U2_[i] / U1_[i]) ** 2)
        F1_[i] = U2_[i]
        F2_[i] = U2_[i] ** 2 / U1_[i] + (g - 1) / g * (U3_[i] - g / 2 * U2_[i] ** 2 / U1_[i])
        # F3_[i] = g * U2_[i] * U3_[i] / U1_[i] - g*(g-1)/2 * U2_[i]**3 / U1_[i]**2
        F3_[i] = g * U2_[i] * U3_[i] / U1_[i] - g * (g - 1) / 2 * (U2_[i] / U1_[i]) ** 2 * U2_[i]

        # Corrector Step
    for i in range(1, len(x) - 1):
        dU1_ = -(F1_[i] - F1_[i - 1]) / dx
        dU2_ = -(F2_[i] - F2_[i - 1]) / dx + 1 / g * r_[i] * T_[i] * (A[i] - A[i - 1]) / dx
        dU3_ = -(F3_[i] - F3_[i - 1]) / dx
        dU1_av = 0.5 * (dU1[i] + dU1_)
        dU2_av = 0.5 * (dU2[i] + dU2_)
        dU3_av = 0.5 * (dU3[i] + dU3_)

        U1[i] = U1[i] + dU1_av * Dt
        U2[i] = U2[i] + dU2_av * Dt
        U3[i] = U3[i] + dU3_av * Dt

        r[i] = U1[i] / A[i]
        V[i] = U2[i] / U1[i]
        T[i] = (g - 1) * (U3[i] / U1[i] - g / 2 * V[i] ** 2)

        # Boundary Conditions
        U2[0] = 2 * U2[1] - U2[2]
        V[0] = U2[0] / U1[0]
        U3[0] = U1[0] * (T[0] / (g - 1) + g / 2 * V[0] ** 2)
        U1[-1] = 2 * U1[-2] - U1[-3]
        U2[-1] = 2 * U2[-2] - U2[-3]
        U3[-1] = 2 * U3[-2] - U3[-3]
        r[-1] = 2 * r[-2] - r[-3]
        V[-1] = 2 * V[-2] - V[-3]
        T[-1] = 2 * T[-2] - T[-3]

    m = r * V * A

    F1 = U2
    F2 = U2 ** 2 / U1 + (g - 1) / g * (U3 - g / 2 * U2 ** 2 / U1)
    F3 = g * U2 * U3 / U1 - g * (g - 1) / 2 * U2 ** 3 / U1 ** 2

    p = r*T
    M = V/(T**2)

    if j == 0 or j == 49 or j == 99 or j == 149 or j == 699:
        mass[j] = m

    r_history.append(r[15])
    T_history.append(T[15])
    p_history.append(p[15])
    M_history.append(M[15])



# Tab. 7.3
print(x.T, A.T, r[:].T, V[:].T, T[:].T, p[:].T, M[:].T, m[:].T)

# Tab. 7.4
print(np.round(np.array((x.T, A.T, r[:].T, r_an[:], np.abs(r[:]-r_an[:])/r[:]*100, M[:].T, Mtot[:].T, np.abs(M[:]-Mtot[:])/M[:]*100)),3))

# Tab. 7.5
# results can be found on Tab. 7.6 if the grid points are changed

# # # Tab. 7.6
# print(f"Density numerical = {r[1399,int(mid)]}, Density analytical = {r_an[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")
# print(f"Temperature numerical = {T[1399,int(mid)]}, Temperature analytical = {T_an[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")
# print(f"Pressure numerical = {p[1399,int(mid)]}, Pressure analytical = {p_an[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")
# print(f"Mach numerical = {M[1399,int(mid)]}, Mach analytical = {Mtot[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")



#
# # Fig. 7.9
plt.figure(figsize=(10, 6))
plt.plot(range(1, len(r_history) + 1), r_history, label="Numerical Solution")
plt.xlabel("Number of Time Steps")
plt.ylabel(r"$\rho / \rho_0$")
plt.title(f"Density at Grid Point {int(mid+1)} Over Time")
plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(range(1, len(T_history) + 1), T_history, label="Numerical Solution")
plt.xlabel("Number of Time Steps")
plt.ylabel(r"$T / T_0$")
plt.title(f"Temperature at Grid Point {int(mid+1)} Over Time")
plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(range(1, len(p_history) + 1), p_history, label="Numerical Solution")
plt.xlabel("Number of Time Steps")
plt.ylabel(r"$p / p_0$")
plt.title(f"Pressure at Grid Point {int(mid+1)} Over Time")
plt.legend()
plt.grid()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(range(1, len(M_history) + 1), M_history, label="Numerical Solution")
plt.xlabel("Number of Time Steps")
plt.ylabel(r"$M$")
plt.title(f"Mach number at Grid Point {int(mid+1)} Over Time")
plt.legend()
plt.grid()
plt.show()
#
# Fig. 7.10
# plt.figure(figsize=(10, 6))
# plt.plot(range(1, len(drdt_av_history) + 1), drdt_av_history, label=r"$|(\frac{d\r}{dt})_{avg}|$")
# plt.plot(range(1, len(dVdt_av_history) + 1), dVdt_av_history, label=r"$|(\frac{dV}{dt})_{avg}|$")
# plt.xlabel("Number of Time Steps")
# plt.ylabel("Residual")
# plt.title(f"Dimensionless time derivatives at Grid Point {int(mid+1)} Over Time")
# plt.legend()
# plt.grid()
# plt.show()
#


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
plt.title("Mass flow distributions")
plt.legend()
plt.grid()
plt.show()
#
# # Fig. 7.12
fig, ax1 = plt.subplots()
ax1.plot(x, r, label="Numerical results")
ax1.plot(x[::3], r_an[::3], 'o', label="Analytical results")
ax1.set_xlabel("Nondimensionless distance through nozzle (x)")
ax1.set_ylabel("Nondimensionless density")
ax2 = ax1.twinx()
ax2.plot(x, M, 'g', label="Numerical results")
ax2.plot(x[::3], Mtot[::3],'o', label="Analytical results")
ax2.set_ylabel("Mach number")
plt.title("Steady-state distributions over nozzle distance")
plt.legend()
plt.grid()
plt.show()

print(mass)