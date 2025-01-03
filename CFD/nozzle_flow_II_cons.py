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
Nt = 5000
x = np.linspace(0,L,Nx) # x/L
dx = L/(Nx-1)
C = 0.5

p0 = 1
pe = 0.93
# pe = 0.9


# initials
rho = 1 - 0.023*x # rho/rho_0
T = 1 - 0.009333*x # T/T_0
V = 0.05+0.11*x # V/a_0
p = (pe-p0)/L *x + 1

# nozzle geometry
limit1 = np.where(x == 1.5)[0][0]
A = np.zeros(len(x))
for i in range(0,limit1):
    A[i] = 1+2.2*(x[i] - 1.5)**2 # A/A*
for i in range(limit1, len(x)):
    A[i] = 1 + 0.2223 * (x[i] - 1.5) ** 2

# Analytical calculations
def Mach_eq(M_ex, pe):
    return pe - (1 + (gamma - 1) / 2 * M_ex ** 2) ** -(gamma / (gamma - 1))


def A_eq(AeA0, M_ex):
    return (AeA0) ** 2 - (1 / M_ex ** 2) * (
                ((2 / (gamma + 1)) * (1 + ((gamma - 1) / 2) * M_ex ** 2)) ** ((gamma + 1) / (gamma - 1)))


mach_ex_guess = 1
A_guess = 0.1
M_ex = fsolve(Mach_eq, mach_ex_guess, args=pe)

AeA0 = fsolve(A_eq, A_guess, args=M_ex)

AA0 = A * AeA0 / 1.5


def Mach_eq2(M_an, AA0):
    return (AA0) ** 2 - (1 / M_an ** 2) * (
                ((2 / (gamma + 1)) * (1 + ((gamma - 1) / 2) * M_an ** 2)) ** ((gamma + 1) / (gamma - 1)))


Mtot = np.zeros(Nx)

for i in range(0, Nx):
    init_guess = 0.079
    M_an = fsolve(Mach_eq2, init_guess, args=AA0[i])

    while M_an < 0.07 or M_an > 0.545:
        init_guess = init_guess + 0.001
        M_an = fsolve(Mach_eq2, init_guess, args=AA0[i])

    Mtot[i] = M_an

p_an = (1 + (gamma - 1) / 2 * Mtot ** 2) ** (-gamma / (gamma - 1))
rho_an = (1 + (gamma - 1) / 2 * Mtot ** 2) ** (-1 / (gamma - 1))
T_an = (1 + (gamma - 1) / 2 * Mtot ** 2) ** -1


U1 = rho * A
U2 = rho * A * V
U3 = rho * (T/(gamma-1) + gamma/2*V**2) * A

F1 = U2
F2 = U2**2/U1 + (gamma-1)/gamma * (U3 - gamma/2 * U2**2/U1)
F3 = gamma * U2*U3/U1 - gamma*(gamma-1)/2 * U2**3/U1**2

J2 = np.zeros(len(x)-1)
dU1 = np.zeros(len(x)-1)
dU2 = np.zeros(len(x)-1)
dU3 = np.zeros(len(x)-1)
Dt = min(C * dx / (T[:]**0.5 + V[:]))
r_ = np.zeros(len(x)-1)
T_ = np.zeros(len(x)-1)

U1_ = np.zeros(len(x)-1)
U2_ = np.zeros(len(x)-1)
U3_ = np.zeros(len(x)-1)
F1_ = np.zeros(len(x)-1)
F2_ = np.zeros(len(x)-1)
F3_ = np.zeros(len(x)-1)
#
# dU1_ = np.zeros(len(x)-1)
# dU2_ = np.zeros(len(x)-1)
# dU3_ = np.zeros(len(x)-1)

for j in range(Nt):
    for i in range (len(x)-1):
        # Predictor Step
        # J2[i] = (gamma-1)/gamma * (U3[i] - gamma/2 * U2[i]**2/U1[i]) * (A[i]-A[i-1])/dx
        J2[i] = 1/gamma * rho[i] * T[i] * (A[i+1]-A[i])/dx
        dU1[i] = -(F1[i+1]-F1[i])/dx
        dU2[i] = -(F2[i+1]-F2[i])/dx + J2[i]
        dU3[i] = -(F3[i+1]-F3[i])/dx
        U1_[i] = U1[i] + dU1[i] * Dt
        U2_[i] = U2[i] + dU2[i] * Dt
        U3_[i] = U3[i] + dU3[i] * Dt
        r_[i] = U1_[i] / A[i]
        T_[i] = (gamma-1) * (U3_[i]/U1_[i] - gamma/2 * (U2_[i]/U1_[i])**2)
        F1_[i] = U2_[i]
        F2_[i] = U2_[i]**2 / U1_[i] + (gamma-1)/gamma * (U3_[i] - gamma/2 * U2_[i]**2 / U1_[i])
        # F3_[i] = gamma * U2_[i] * U3_[i] / U1_[i] - gamma*(gamma-1)/2 * U2_[i]**3 / U1_[i]**2
        F3_[i] = gamma * U2_[i] * U3_[i] / U1_[i] - gamma*(gamma-1)/2 * (U2_[i]/U1_[i])**2 * U2_[i]

    # print(dU1_[[15,14]])
    # print(dU2_[[15,14]])
    # print(dU3[[15,14]])

        # Corrector Step
    for i in range (1,len(x)-1):
        dU1_ = -(F1_[i]-F1_[i-1])/dx
        dU2_ = -(F2_[i]-F2_[i-1])/dx + 1/gamma * r_[i] * T_[i] * (A[i]-A[i-1])/dx
        dU3_ = -(F3_[i]-F3_[i-1])/dx
        dU1_av = 0.5 * (dU1[i] + dU1_)
        dU2_av = 0.5 * (dU2[i] + dU2_)
        dU3_av = 0.5 * (dU3[i] + dU3_)

        U1[i] = U1[i] + dU1_av * Dt
        U2[i] = U2[i] + dU2_av * Dt
        U3[i] = U3[i] + dU3_av * Dt

        rho[i] = U1[i]/A[i]
        V[i] = U2[i]/U1[i]
        T[i] = (gamma-1) * (U3[i]/U1[i] - gamma/2 * V[i]**2)

        # Boundary Conditions
        U2[0] = 2 * U2[1] - U2[2]
        V[0] = U2[0]/U1[0]
        U3[0] = U1[0] * (T[0]/(gamma-1) + gamma/2 * V[0]**2)
        U1[-1] = 2 * U1[-2] - U1[-3]
        U2[-1] = 2 * U2[-2] - U2[-3]
        V[-1] = 2 * V[-2] - V[-3]
        rho[-1] = 2 * rho[-2] - rho[-3]
        T[-1] = pe/rho[-1]
        U3[-1] = U1[-1] * (T[-1]/(gamma-1) + gamma/2 * V[-1]**2)


        m = rho * V * A
    # print(rho[15])
    # print(T[15])
    # print(V[15])

        # if i ==15:
        #     print(dU1_)
        #     print(dU2_av)
        #     print(dU3_av)

    F1 = U2
    F2 = U2 ** 2 / U1 + (gamma - 1) / gamma * (U3 - gamma / 2 * U2 ** 2 / U1)
    F3 = gamma * U2 * U3 / U1 - gamma * (gamma - 1) / 2 * U2 ** 3 / U1 ** 2

    p = rho*T
    M = V/(T**2)

    if j == 0 or j == 499 or j == Nt-1:
        mass[j] = m
        pressure[j] = p
    if j == 999:
        pressure[j] = p
    if j == 399 or j == 799 or j == 2199:
        pressure[j] = p


print(rho)
print(V)
print(T)
print(m)


pressure[0] = (pe-p0)/L *x + 1

########## RESULTS ###########

# Tab. 7.7
print(x.T, A.T, rho[:].T, V[:].T, T[:].T, p[:].T, M[:].T, m[:   ].T)

# Tab. 7.8
print(np.round(np.array((x.T, A.T, rho[:].T, rho_an[:], np.abs(rho[:]-rho_an[:])/rho[:]*100, M[:].T, Mtot[:].T, np.abs(M[:]-Mtot[:])/M[:]*100)),3))

# Tab. 7.5
# results can be found on Tab. 7.6 if the grid points are changed

# # Tab. 7.6
# print(f"Density numerical = {rho[1399,int(mid)]}, Density analytical = {rho_an[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")
# print(f"Temperature numerical = {T[1399,int(mid)]}, Temperature analytical = {T_an[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")
# print(f"Pressure numerical = {p[1399,int(mid)]}, Pressure analytical = {p_an[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")
# print(f"Mach numerical = {M[1399,int(mid)]}, Mach analytical = {Mtot[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")



#
# # Fig. 7.9
# plt.figure(figsize=(10, 6))
# plt.plot(range(1, len(rho_history) + 1), rho_history, label="Numerical Solution")
# plt.xlabel("Number of Time Steps")
# plt.ylabel(rho"$\rho / \rho_0$")
# plt.title(f"Density at Grid Point {int(mid+1)} Over Time")
# plt.legend()
# plt.grid()
# plt.show()
#
# plt.figure(figsize=(10, 6))
# plt.plot(range(1, len(T_history) + 1), T_history, label="Numerical Solution")
# plt.xlabel("Number of Time Steps")
# plt.ylabel(rho"$T / T_0$")
# plt.title(f"Temperature at Grid Point {int(mid+1)} Over Time")
# plt.legend()
# plt.grid()
# plt.show()
#
# plt.figure(figsize=(10, 6))
# plt.plot(range(1, len(p_history) + 1), p_history, label="Numerical Solution")
# plt.xlabel("Number of Time Steps")
# plt.ylabel(rho"$p / p_0$")
# plt.title(f"Pressure at Grid Point {int(mid+1)} Over Time")
# plt.legend()
# plt.grid()
# plt.show()
#
# plt.figure(figsize=(10, 6))
# plt.plot(range(1, len(M_history) + 1), M_history, label="Numerical Solution")
# plt.xlabel("Number of Time Steps")
# plt.ylabel(rho"$M$")
# plt.title(f"Mach number at Grid Point {int(mid+1)} Over Time")
# plt.legend()
# plt.grid()
# plt.show()
#
# Fig. 7.10
# plt.figure(figsize=(10, 6))
# plt.plot(range(1, len(drhodt_av_history) + 1), drhodt_av_history, label=rho"$|(\frac{d\rho}{dt})_{avg}|$")
# plt.plot(range(1, len(dVdt_av_history) + 1), dVdt_av_history, label=rho"$|(\frac{dV}{dt})_{avg}|$")
# plt.xlabel("Number of Time Steps")
# plt.ylabel("Residual")
# plt.title(f"Dimensionless time derivatives at Grid Point {int(mid+1)} Over Time")
# plt.legend()
# plt.grid()
# plt.show()
#


# # Fig. 7.16
plt.figure(figsize=(10, 6))
plt.plot(x,mass[0], label=r"$0\Delta t$")
plt.plot(x,mass[499], label=r"$500\Delta t$")
plt.plot(x,mass[Nt-1], label=r"$5000\Delta t$")
plt.xlabel("Nondimensionless distance through nozzle (x)")
plt.ylabel("Nondimensionless mass flow")
plt.title("Mass flow distributions")
plt.legend()
plt.grid()
plt.show()
#
# # Fig. 7.17
plt.figure(figsize=(10, 6))
plt.plot(x,pressure[0], label=r"$0\Delta t$")
plt.plot(x,pressure[499], label=r"$500\Delta t$")
plt.plot(x,pressure[999], label=r"$1000\Delta t$")
plt.plot(x,pressure[Nt-1], label=r"$5000\Delta t$")
plt.xlabel("Nondimensionless distance through nozzle (x)")
plt.ylabel("Nondimensionless mass flow")
plt.title("Mass flow distributions")
plt.legend()
plt.grid()
plt.show()

# Fig. 7.18
############### CHANGE pe to pe = 0.9 for this plot
# plt.figure(figsize=(10, 6))
# plt.plot(x,pressure[0], label=r"$0\Delta t$")
# plt.plot(x,pressure[399], label=r"$400\Delta t$")
# plt.plot(x,pressure[799], label=r"$800\Delta t$")
# plt.plot(x,pressure[2199], label=r"$1200\Delta t$")
# plt.xlabel("Nondimensionless distance through nozzle (x)")
# plt.ylabel("Nondimensionless mass flow")
# plt.title("Mass flow distributions")
# plt.legend()
# plt.grid()
# plt.show()
