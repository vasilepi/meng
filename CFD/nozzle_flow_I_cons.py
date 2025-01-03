import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import fsolve


C = 1.1  # Courant
g = 1.4
L = 3  # m
Nx = 31
dx = L / (Nx-1)
Nt = 1400
x = np.linspace(0, L, Nx)
A = 1 + 2.2 * (x - 1.5) ** 2
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

    for i in range(1,Nx-1):
        # J2_est[t+1,i] = 1 / gamma * rho_est[t + 1, i] * T_est[t + 1, i] * (A[i+1] - A[i]) / dx
        J2_est[t + 1, i] = (gamma-1) / gamma * (U3_est[t + 1, i] * 0.5*gamma*U2_est[t + 1, i]**2/U1_est[t+1,i]) * (np.log(A[i]) - np.log(A[i-1])) / dx
    for i in range(1,Nx-1):
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



