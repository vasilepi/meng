from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import minimize

P = 5000 # N
A1 = 1.8  # Area bar 1 [mm^2]
A2 = 3.1  # Area bar 2 [mm^2]
l1 = 500  # Length bar 1 [mm]
l2 = 1000 # Length bar 2 [mm]
l = 1505  # Distance between walls [mm]
E = 210e3 # Young's modulus [MPa]
e = 1e10   # Penalty stiffness [N/mm]

# Stiffnesses
k1 = E * A1 / l1
k2 = E * A2 / l2

K1 = k1 * np.array([
    [1, -1],
    [-1, 1]
])

K2 = k2 * np.array([
    [1, -1],
    [-1, 1]
])


K = np.zeros((4,4))
K[0:2, 0:2] += K1
K[2:4, 2:4] += K2



# F = np.array([0., P, 0., 0.])


# BCs
fix = [0, 4]
free = [i for i in range(0,3) if i not in fix]
K = K[np.ix_(free, free)]
# F = F[free]


Farr = np.linspace(0,P,100)

u1=[]
u2=[]
N = []
gap = []

g0 = l - l1 - l2

epsilon = np.eye(K.shape[0]) * e
# print(epsilon)
for f in Farr:
    F = np.array([0,f,0,0])
    F = F[free]
    # print(F)
    u = np.linalg.solve(K,F)

    g = g0 - u
    # print(g)
    if g[0] > 0:
        u1.append(u[0])
        u2.append(u[1])
        N.append(0)
        gap.append(g[0])
    # print(u)

    else:
        Knew = np.linalg.inv(epsilon)*K + np.eye(K.shape[0])
        # print(F)
        Fnew = np.linalg.inv(epsilon) @ F
        # print("K",Knew)
        # print("F",Fnew)
        u = np.linalg.solve(Knew, Fnew)
        # print(u)
        
        # print(gnew)
        u1.append(u[0] + u1[-1])
        u2.append(u2[-1]-u[0])
        gnew1 = u1[-1] - g0
        gnew2 = g0 - u2[-1]
        gnew = np.array([gnew1, gnew2])
        FN = -epsilon @ gnew
        # pri
        N.append(FN[0])
        gap.append(gap[-1] + gnew[0])
    print("GAP", gap[-1])

s1 = np.array(u1) * E/l1
s2 = np.array(u2) * E/l2


sf = 10000

print("DISP B1",u1[-1])
print("STRESS B1", s1[-1])
print("DISP B2", u2[-1])
print("STRESS B2", s2[-1])
print("CONT FORCE", N[-1])


plt.figure()
plt.plot(u1, Farr)
plt.xlabel("Displacement of Node 1 (mm)")
plt.ylabel("External Force (N)")
plt.grid()
plt.show()

plt.figure()
plt.plot(Farr, s1, label = "Beam 1")
plt.plot(Farr, s2, label = "Beam 2")
plt.ylabel("Principal Stress (MPa)")
plt.xlabel("External Force (N)")
plt.legend()
plt.grid()
plt.show()

plt.figure()
plt.plot(Farr, s1, label = "Beam 1")
plt.plot(Farr, s2*sf, label = "Beam 2 - Scaled")
plt.ylabel("Principal Stress (MPa)")
plt.xlabel("External Force (N)")
plt.legend()
plt.grid()
plt.show()

plt.figure()
plt.plot(Farr, N)
plt.ylabel("Contact Force (N)")
plt.xlabel("External Force (N)")
plt.grid()
plt.show()


plt.figure()
plt.plot(Farr, u1, label = "Beam 1")
plt.plot(Farr, u2, label = "Beam 1")
plt.plot(Farr, gap, label = "Gap")
plt.ylabel("Displacement (mm)")
plt.xlabel("External Force (N)")
plt.legend()
plt.grid()
plt.show()