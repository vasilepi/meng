import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve




def conservation_check(rho, V, T, A, t, dx):
    mass = np.sum(rho[t, :] * A) * dx
    momentum = np.sum(rho[t, :] * V[t, :] * A) * dx
    energy = np.sum(rho[t, :] * T[t, :] * A) * dx
    return mass, momentum, energy


gamma = 1.4
L = 3
Nx = 61
Nt = 1600
x = np.linspace(0,L,Nx) # x/L
dx = L/(Nx-1)
C = 0.5
Cx = 0.2
p0 = 1
pe = 0.6784

rho = np.zeros(len(x))
T = np.zeros(len(x))
p = np.zeros(len(x))
# initials

limit1 = np.where(x >= 0.5)[0][0]
limit2 = np.where(x >= 1.5)[0][0]
limit3 = np.where(x >= 2.1)[0][0]
for i in range(0,limit1):
    rho[i] = 1
    T[i] = 1
for i in range(limit1,limit2):
    rho[i] = 1-0.366*(x[i]-0.5)
    T[i] = 1 - 0.167 * (x[i] - 0.5)
for i in range(limit2,limit3):
    rho[i] = 0.634-0.702*(x[i]-1.5)
    T[i] = 0.833 - 0.4908*(x[i]-1.5)
for i in range(limit3,Nx):
    rho[i] = 0.5892+0.10228*(x[i]-2.1)
    T[i] = 0.93968+0.0622*(x[i]-2.1)

# nozzle geometry
A_ = 1 + 2.2 * (x - 1.5) ** 2
A = A_/min(A_)
Ae = A[-1]
V = 0.59 / (rho * A)
p[:] = rho[:] * T[:]


# Analytical calculations 

def Me_solve(Me, pe):
    return pe * Ae - (1 + (gamma - 1) / 2 * Me ** 2) ** -(gamma / (gamma - 1)) * (1 / Me) * (
                ((2 / (gamma + 1)) * (1 + ((gamma - 1) / 2) * Me ** 2)) ** ((gamma + 1) / (2 * (gamma - 1))))


guess = 0.1
Me = fsolve(Me_solve, guess, args=pe)

pep0e = (1 + (gamma - 1) / 2 * Me ** 2) ** -(gamma / (gamma - 1))
p02p01 = pe / pep0e

TeT0 = (1 + (gamma - 1) / 2 * Me ** 2) ** -1
rhoerho0 = (1 + (gamma - 1) / 2 * Me ** 2) ** -(1 / (gamma - 1)) * p02p01
me = rhoerho0 * Me * Ae


Mt = 1
pt = (1 + (gamma - 1) / 2 * Mt ** 2) ** -(gamma / (gamma - 1))
rhot = (1 + (gamma - 1) / 2 * Mt ** 2) ** -(1 / (gamma - 1))
Tt = (1 + (gamma - 1) / 2 * Mt ** 2) ** -1

def M1_solve(M1, p02p01):
    return p02p01 - ((gamma + 1) * M1 ** 2 / ((gamma - 1) * M1 ** 2 + 2)) ** (gamma / (gamma - 1)) * (
                (gamma + 1) / (2 * gamma * M1 ** 2 - gamma + 1)) ** (1 / (gamma - 1))


M1_guess = 2
M1 = fsolve(M1_solve, M1_guess, args=p02p01)


def Astar_solve(AAdot, M):
    return (AAdot) ** 2 - (1 / M ** 2) * (
                ((2 / (gamma + 1)) * (1 + ((gamma - 1) / 2) * M ** 2)) ** ((gamma + 1) / (gamma - 1)))


A1A1star = fsolve(Astar_solve, 1, args=M1)

p2p1 = 1 + 2 * gamma / (gamma + 1) * (M1 ** 2 - 1)
M2 = ((1 + (gamma - 1) / 2 * M1 ** 2) / (gamma * M1 ** 2 - (gamma - 1) / 2)) ** 0.5




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
p_est = np.zeros(len(x))
V_est = np.zeros(len(x))
U1_est = np.zeros(len(x))
U2_est = np.zeros(len(x))
U3_est = np.zeros(len(x))
F1_est = np.zeros(len(x)-1)
F2_est = np.zeros(len(x)-1)
F3_est = np.zeros(len(x)-1)

S1 = np.zeros(Nx)
S2 = np.zeros(Nx)
S3 = np.zeros(Nx)
S1_est = np.zeros(Nx)
S2_est = np.zeros(Nx)
S3_est = np.zeros(Nx)


mid = (Nx-1)/2
mass = {}
pressure = {}
for j in range(Nt):
    for i in range(1, len(x) - 1):
        S1[i] = Cx * abs(p[i + 1] - 2 * p[i] + p[i - 1]) / (p[i + 1] + 2 * p[i] + p[i - 1]) * (
                    U1[i + 1] - 2 * U1[i] + U1[i - 1])
        S2[i] = Cx * abs(p[i + 1] - 2 * p[i] + p[i - 1]) / (p[i + 1] + 2 * p[i] + p[i - 1]) * (
                U2[i + 1] - 2 * U2[i] + U2[i - 1])
        S3[i] = Cx * abs(p[i + 1] - 2 * p[i] + p[i - 1]) / (p[i + 1] + 2 * p[i] + p[i - 1]) * (
                U3[i + 1] - 2 * U3[i] + U3[i - 1])

        ########### For Fig. 7.23 ############
        # S1[i] = 0
        # S2[i] = 0
        # S3[i] = 0
    for i in range (len(x)-1):
        # Predictor Step
        J2[i] = 1/gamma * rho[i] * T[i] * (A[i+1]-A[i])/dx
        dU1dt[i] = -(F1[i+1]-F1[i])/dx
        dU2dt[i] = -(F2[i+1]-F2[i])/dx + J2[i]
        dU3dt[i] = -(F3[i+1]-F3[i])/dx
        U1_est[i] = U1[i] + dU1dt[i] * dt + S1[i]
        U2_est[i] = U2[i] + dU2dt[i] * dt + S2[i]
        U3_est[i] = U3[i] + dU3dt[i] * dt + S3[i]

        rho_est[i] = U1_est[i] / A[i]
        T_est[i] = (gamma-1) * (U3_est[i]/U1_est[i] - gamma/2 * (U2_est[i]/U1_est[i])**2)
        p_est[i] = rho_est[i] * T_est[i]
        F1_est[i] = U2_est[i]
        F2_est[i] = U2_est[i]**2 / U1_est[i] + (gamma-1)/gamma * (U3_est[i] - gamma/2 * U2_est[i]**2 / U1_est[i])
        # F3_est[i] = gamma * U2_est[i] * U3_est[i] / U1_est[i] - gamma*(gamma-1)/2 * U2_est[i]**3 / U1_est[i]**2
        F3_est[i] = gamma * U2_est[i] * U3_est[i] / U1_est[i] - gamma*(gamma-1)/2 * (U2_est[i]/U1_est[i])**2 * U2_est[i]
    U1_est[-1] = 2 * U1_est[-2] - U1_est[-3]
    U2_est[-1] = 2 * U2_est[-2] - U2_est[-3]
    V_est = U2_est/U1_est
    U3_est[-1] = pe * A[-1] / (gamma - 1) + 0.5 * gamma * U2_est[-1] * V_est[-1]
    p_est[0] = 1
    p_est[-1] = pe
    for i in range (1, len(x)-1):
        S1_est[i] = Cx * abs(p_est[i + 1] - 2 * p_est[i] + p_est[i - 1]) / (p_est[i + 1] + 2 * p_est[i] + p_est[i - 1]) * (
                U1_est[i + 1] - 2 * U1_est[i] + U1_est[i - 1])
        S2_est[i] = Cx * abs(p_est[i + 1] - 2 * p_est[i] + p_est[i - 1]) / (
                    p_est[i + 1] + 2 * p_est[i] + p_est[i - 1]) * (
                            U2_est[i + 1] - 2 * U2_est[i] + U2_est[i - 1])
        S3_est[i] = Cx * abs(p_est[i + 1] - 2 * p_est[i] + p_est[i - 1]) / (
                    p_est[i + 1] + 2 * p_est[i] + p_est[i - 1]) * (
                            U3_est[i + 1] - 2 * U3_est[i] + U3_est[i - 1])
        ########### For Fig. 7.23 ############
        # S1_est[i] = 0
        # S2_est[i] = 0
        # S3_est[i] = 0
        # Corrector Step
    for i in range (1,len(x)-1):
        dU1dt_est = -(F1_est[i]-F1_est[i-1])/dx
        dU2dt_est = -(F2_est[i]-F2_est[i-1])/dx + 1/gamma * rho_est[i] * T_est[i] * (A[i]-A[i-1])/dx
        dU3dt_est = -(F3_est[i]-F3_est[i-1])/dx
        dU1dt_av = 0.5 * (dU1dt[i] + dU1dt_est)
        dU2dt_av = 0.5 * (dU2dt[i] + dU2dt_est)
        dU3dt_av = 0.5 * (dU3dt[i] + dU3dt_est)

        U1[i] = U1[i] + dU1dt_av * dt + S1_est[i]
        U2[i] = U2[i] + dU2dt_av * dt + S2_est[i]
        U3[i] = U3[i] + dU3dt_av * dt + S3_est[i]

        rho[i] = U1[i]/A[i]
        V[i] = U2[i]/U1[i]
        T[i] = (gamma-1) * (U3[i]/U1[i] - gamma/2 * V[i]**2)
        p[i] = rho[i]*T[i]


        # Boundary Conditions
        U1[0] = A[0]
        U2[0] = 2 * U2[1] - U2[2]
        V[0] = U2[0] / U1[0]
        U3[0] = U1[0] * (T[0] / (gamma - 1) + gamma / 2 * V[0] ** 2)
        rho[0] = 1
        T[0] = 1
        p[0] = 1

        U1[-1] = 2 * U1[-2] - U1[-3]
        U2[-1] = 2 * U2[-2] - U2[-3]
        V[-1] = U2[-1] / U1[-1]
        U3[-1] = pe*A[-1]/(gamma-1)+0.5*gamma*U2[-1]*V[-1]
        rho[-1] = U1[-1]/A[-1]
        T[-1] = pe/rho[-1]
        p[-1] = pe


    F1 = U2
    F2 = U2 ** 2 / U1 + (gamma - 1) / gamma * (U3 - gamma / 2 * U2 ** 2 / U1)
    F3 = gamma * U2 * U3 / U1 - gamma * (gamma - 1) / 2 * U2 ** 3 / U1 ** 2

    p = rho*T
    M = V/(T**0.5)
    mass_flow = rho * V * A



    # if j == 0 or j == 499 or j == Nt-1:
    #     mass[j] = mass_flow
    #     pressure[j] = p
    # if j == 999:
    #     pressure[j] = p
    # if j == 399 or j == 799 or j == 2199:
    #     pressure[j] = p




########## RESULTS ###########

# Fig. 7.23 - 7.24 (For 7.23 replace viscosity with commented values)
plt.figure(figsize=(10, 6))
plt.plot(x,p, label="$Numerical$")
plt.xlabel("Nondimensionless distance through nozzle (x)")
plt.ylabel("Nondimensionless density")
plt.legend()
plt.grid()
plt.show()

# Fig. 7.25
plt.figure(figsize=(10, 6))
plt.plot(x,M, label=r"$0\Delta t$")
plt.xlabel("Nondimensionless distance through nozzle (x)")
plt.ylabel("Mach")
plt.legend()
plt.grid()
plt.show()

# Fig. 7.26
plt.figure(figsize=(10, 6))
plt.plot(x,mass_flow, label=r"$0\Delta t$")
plt.xlabel("Nondimensionless distance through nozzle (x)")
plt.ylabel("Nondimensionless mass flow")
# plt.title("Mass flow distributions")
plt.legend()
plt.grid()
plt.show()


# Tab. 7.13
print(x.T, A.T, rho[:].T, V[:].T, T[:].T, p[:].T, M[:].T, mass_flow[:].T)

# Tab. 7.8
# print(np.round(np.array((x.T, A.T, rho[:].T, rho_an[:], np.abs(rho[:]-rho_an[:])/rho[:]*100, M[:].T, M_an[:].T, np.abs(M[:]-M_an[:])/M[:]*100)),3))

# Tab. 7.5
# results can be found on Tab. 7.6 if the grid points are changed

# # Tab. 7.14
print("------ TAB 7.14 --------")
print(f"Density numerical = {rho[int(mid)]}, Density analytical = {rhot} for C = {C} at throat")
print(f"Temperature numerical = {T[int(mid)]}, Temperature analytical = {Tt} for C = {C} at throat")
print(f"Pressure numerical = {p[int(mid)]}, Pressure analytical = {pt} for C = {C} at throat")
print(f"Mach numerical = {M[int(mid)]}, Mach analytical = {Mt} for C = {C} at throat")
print(f"Mass numerical = {mass_flow[int(mid)]}, Mass analytical = {me} for C = {C} at throat")


# # Tab. 7.15
print("------ TAB 7.15 --------")
print(f"Density numerical = {rho[-1]}, Density analytical = {rhoerho0} for C = {C} at exit")
print(f"Temperature numerical = {T[-1]}, Temperature analytical = {TeT0} for C = {C} at exit")
print(f"Pressure numerical = {p[-1]}, Pressure analytical = {pe} for C = {C} at exit")
print(f"Mach numerical = {M[-1]}, Mach analytical = {Me} for C = {C} at exit")
print(f"Mass numerical = {mass_flow[-1]}, Mass analytical = {me} for C = {C} at exit")

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
# plt.figure(figsize=(10, 6))
# plt.plot(x,mass[0], label=r"$0\Delta t$")
# plt.plot(x,mass[499], label=r"$500\Delta t$")
# plt.plot(x,mass[Nt-1], label=r"$5000\Delta t$")
# plt.xlabel("Nondimensionless distance through nozzle (x)")
# plt.ylabel("Nondimensionless mass flow")
# plt.title("Mass flow distributions")
# plt.legend()
# plt.grid()
# plt.show()
#
# # Fig. 7.17
# plt.figure(figsize=(10, 6))
# plt.plot(x,pressure[0], label=r"$0\Delta t$")
# plt.plot(x,pressure[499], label=r"$500\Delta t$")
# plt.plot(x,pressure[999], label=r"$1000\Delta t$")
# plt.plot(x,pressure[Nt-1], label=r"$5000\Delta t$")
# plt.xlabel("Nondimensionless distance through nozzle (x)")
# plt.ylabel("Nondimensionless mass flow")
# plt.title("Mass flow distributions")
# plt.legend()
# plt.grid()
# plt.show()

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

