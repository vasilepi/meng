import numpy as np
import matplotlib as plt


gamma = 1.4
Nx = 31
x = np.linspace(0,3,Nx) # x/L
dx = 3/(Nx-1)
tss = 500
C = 0.5

# initials
rho_i = 1 - 0.3146*x # rho/rho_0
T_i = 1 - 0.2314*x # T/T_0
V_i = (0.1 + 1.09*x)*T_i**0.5 # V/a_0

# nozzle geometry
A = 1+2.2*(x - 1.5)**2 # A/A*

dt = np.min(C*dx/(T_i**0.5+V_i))
Nt = tss/dt + 1
t = np.linspace(0, tss, int(Nt))


rho = np.zeros((len(t),len(x)))
V = np.zeros((len(t),len(x)))
T = np.zeros((len(t),len(x)))
drhodt = np.zeros((len(t),len(x)))
dVdt = np.zeros((len(t),len(x)))
dTdt = np.zeros((len(t),len(x)))

rho[0,:] = rho_i
V[0,:] = V_i
T[0,:] = T_i



rho_est = np.zeros((len(t),len(x)))
V_est = np.zeros((len(t),len(x)))
T_est = np.zeros((len(t),len(x)))
drhodt_est = np.zeros((len(t),len(x)))
dVdt_est = np.zeros((len(t),len(x)))
dTdt_est = np.zeros((len(t),len(x)))


for t in range(1):
    for i in range(1,Nx-2):
        # derivatives
        drhodt[t, i] = - rho[t,i] * (V[t,i+1]-V[t,i])/dx - rho[t,i]*V[t,i]*(np.log(A[i+1])-np.log(A[i]))/dx - V[t,i]*(rho[t,i+1]-rho[t,i])/dx
        dVdt[t,i] = -V[t,i]*(V[t,i+1]-V[t,i])/dx - 1/gamma *((T[t,i+1]-T[t,i])/dx + T[t,i]/rho[t,i]*(rho[t,i+1]-rho[t,i])/dx)
        dTdt[t, i] = -V[t,i]*(T[t,i+1]-T[t,i])/dx - (gamma-1)*T[t,i]*((V[t,i+1]-V[t,i])/dx + V[t,i]*(np.log(A[i+1]) - np.log(A[i]))/dx)

        # estimators
        rho_est[t+1,i] = rho[t, i] + drhodt[t,i]*dt
        V_est[t + 1, i] = V[t, i] + dVdt[t, i] * dt
        T_est[t + 1, i] = T[t, i] + dTdt[t, i] * dt

        # der estimators
        drhodt_est[t+1, i] = - rho_est[t+1, i] * (V_est[t+1, i] - V_est[t+1, i-1]) / dx - rho_est[t+1, i] * V_est[t+1, i] * (np.log(A[i]) - np.log(A[i-1])) / dx - V_est[t+1, i] * (rho_est[t+1, i] - rho_est[t+1, i-1]) / dx
        dVdt_est[t+1, i] = - V_est[t + 1, i] * (V_est[t + 1, i] - V_est[t + 1, i-1]) / dx - 1/gamma* ((T_est[t + 1, i]-T_est[t+1,i-1])/dx + T_est[t + 1, i]/rho_est[t+1,i] * (rho_est[t + 1, i] - rho_est[t + 1, i-1]) / dx)
        dTdt_est[t+1, i] = - V_est[t + 1, i] * (T_est[t + 1, i] - T_est[t + 1, i-1]) / dx - (gamma-1) * T_est[t + 1, i] * ((V_est[t + 1, i] - V_est[t + 1, i-1]) / dx + V_est[t+1,i]*(np.log(A[i]) - np.log(A[i-1])) / dx)

print(drhodt[0,15], dVdt[0,15], dTdt[0,15], rho_est[1,15], V_est[1,15], T_est[1,15], drhodt_est[1,15], dVdt_est[1,15], dTdt_est[1,15])
#     for i in range(1,len(x)):
#         drhodt[0,i] = 1










