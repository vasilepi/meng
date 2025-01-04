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
def M_solve(M_ex, pe):
    return pe - (1 + (gamma - 1) / 2 * M_ex ** 2) ** -(gamma / (gamma - 1))

def A_solve(AeA0, M_ex):
    return (AeA0) ** 2 - (1 / M_ex ** 2) * (
        ((2 / (gamma + 1)) * (1 + ((gamma - 1) / 2) * M_ex ** 2)) ** ((gamma + 1) / (gamma - 1))
    )

mach_ex_guess = 1
A_guess = 0.1
M_ex = fsolve(M_solve, mach_ex_guess, args=pe)[0]  # Extract the scalar
AeA0 = fsolve(A_solve, A_guess, args=M_ex)[0]  # Extract the scalar
AA0 = A * AeA0 / 1.5

def M_solvesolve(M_, AA0):
    return (AA0) ** 2 - (1 / M_ ** 2) * (
        ((2 / (gamma + 1)) * (1 + ((gamma - 1) / 2) * M_ ** 2)) ** ((gamma + 1) / (gamma - 1))
    )

M_an = np.zeros(Nx)

for i in range(0, Nx):
    init_guess = 0.079
    M_ = fsolve(M_solvesolve, init_guess, args=AA0[i])[0]  # Extract the scalar

    while M_ < 0.07 or M_ > 0.545:
        init_guess = init_guess + 0.001
        M_ = fsolve(M_solvesolve, init_guess, args=AA0[i])[0]  # Extract the scalar

    M_an[i] = M_

# Calculate analytical results
p_an = (1 + (gamma - 1) / 2 * M_an ** 2) ** (-gamma / (gamma - 1))
rho_an = (1 + (gamma - 1) / 2 * M_an ** 2) ** (-1 / (gamma - 1))
T_an = (1 + (gamma - 1) / 2 * M_an ** 2) ** -1

U1 = rho * A
U2 = rho * A * V
U3 = rho * (T/(gamma-1) + gamma/2*V**2) * A

F1 = U2
F2 = U2**2/U1 + (gamma-1)/gamma * (U3 - gamma/2 * U2**2/U1)
F3 = gamma * U2*U3/U1 - gamma*(gamma-1)/2 * U2**3/U1**2

J2 = np.zeros(len(x)-1)
dU1dt = np.zeros(len(x)-1)
dU2dt = np.zeros(len(x)-1)
dU3dt = np.zeros(len(x)-1)
dt = min(C * dx / (T[:]**0.5 + V[:]))
rho_est = np.zeros(len(x)-1)
T_est = np.zeros(len(x)-1)

U1_est = np.zeros(len(x)-1)
U2_est = np.zeros(len(x)-1)
U3_est = np.zeros(len(x)-1)
F1_est = np.zeros(len(x)-1)
F2_est = np.zeros(len(x)-1)
F3_est = np.zeros(len(x)-1)

mass = {}
pressure = {}
for j in range(Nt):
    for i in range (len(x)-1):
        # Predictor Step
        # J2[i] = (gamma-1)/gamma * (U3[i] - gamma/2 * U2[i]**2/U1[i]) * (A[i]-A[i-1])/dx
        J2[i] = 1/gamma * rho[i] * T[i] * (A[i+1]-A[i])/dx
        dU1dt[i] = -(F1[i+1]-F1[i])/dx
        dU2dt[i] = -(F2[i+1]-F2[i])/dx + J2[i]
        dU3dt[i] = -(F3[i+1]-F3[i])/dx
        U1_est[i] = U1[i] + dU1dt[i] * dt
        U2_est[i] = U2[i] + dU2dt[i] * dt
        U3_est[i] = U3[i] + dU3dt[i] * dt
        rho_est[i] = U1_est[i] / A[i]
        T_est[i] = (gamma-1) * (U3_est[i]/U1_est[i] - gamma/2 * (U2_est[i]/U1_est[i])**2)
        F1_est[i] = U2_est[i]
        F2_est[i] = U2_est[i]**2 / U1_est[i] + (gamma-1)/gamma * (U3_est[i] - gamma/2 * U2_est[i]**2 / U1_est[i])
        # F3_est[i] = gamma * U2_est[i] * U3_est[i] / U1_est[i] - gamma*(gamma-1)/2 * U2_est[i]**3 / U1_est[i]**2
        F3_est[i] = gamma * U2_est[i] * U3_est[i] / U1_est[i] - gamma*(gamma-1)/2 * (U2_est[i]/U1_est[i])**2 * U2_est[i]



        # Corrector Step
    for i in range (1,len(x)-1):
        dU1dt_est = -(F1_est[i]-F1_est[i-1])/dx
        dU2dt_est = -(F2_est[i]-F2_est[i-1])/dx + 1/gamma * rho_est[i] * T_est[i] * (A[i]-A[i-1])/dx
        dU3dt_est = -(F3_est[i]-F3_est[i-1])/dx
        dU1dt_av = 0.5 * (dU1dt[i] + dU1dt_est)
        dU2dt_av = 0.5 * (dU2dt[i] + dU2dt_est)
        dU3dt_av = 0.5 * (dU3dt[i] + dU3dt_est)

        U1[i] = U1[i] + dU1dt_av * dt
        U2[i] = U2[i] + dU2dt_av * dt
        U3[i] = U3[i] + dU3dt_av * dt

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

    mass_flow = rho * V * A



    F1 = U2
    F2 = U2 ** 2 / U1 + (gamma - 1) / gamma * (U3 - gamma / 2 * U2 ** 2 / U1)
    F3 = gamma * U2 * U3 / U1 - gamma * (gamma - 1) / 2 * U2 ** 3 / U1 ** 2

    p = rho*T
    M = V/(T**2)

    if j == 0 or j == 499 or j == Nt-1:
        mass[j] = mass_flow
        pressure[j] = p
    if j == 999:
        pressure[j] = p
    if j == 399 or j == 799 or j == 2199:
        pressure[j] = p


pressure[0] = (pe-p0)/L *x + 1

########## RESULTS ###########

# Tab. 7.7
print(x.T, A.T, rho[:].T, V[:].T, T[:].T, p[:].T, M[:].T, mass_flow[:   ].T)

# Tab. 7.8
print(np.round(np.array((x.T, A.T, rho[:].T, rho_an[:], np.abs(rho[:]-rho_an[:])/rho[:]*100, M[:].T, M_an[:].T, np.abs(M[:]-M_an[:])/M[:]*100)),3))

# Tab. 7.5
# results can be found on Tab. 7.6 if the grid points are changed

# # Tab. 7.6
# print(f"Density numerical = {rho[1399,int(mid)]}, Density analytical = {rho_an[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")
# print(f"Temperature numerical = {T[1399,int(mid)]}, Temperature analytical = {T_an[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")
# print(f"Pressure numerical = {p[1399,int(mid)]}, Pressure analytical = {p_an[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")
# print(f"Mach numerical = {M[1399,int(mid)]}, Mach analytical = {M_an[int(mid)]} for C = {C} at GRID POINT {int(mid+1)}")



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
