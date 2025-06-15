import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# parameters
E = 210000  # MPa
l = 1505  # mm
l1 = 500       # mm
l2 = 1000      # mm
A1 = 1.8         # cm²
A2 = 3.1         # cm²

k1 = E * A1 / l1
k2 = E * A2 / l2

# stiffness matrices
k_elem_1 = k1 * np.array([[1, -1], [-1, 1]])
k_elem_2 = k2 * np.array([[1, -1], [-1, 1]])


K_global = np.zeros((4, 4))
K_global[:2, :2] += k_elem_1
K_global[2:, 2:] += k_elem_2

# sim
initial_gap = l - l1 - l2
penalty = 1e10
steps = 100
max_load = 8000
f_array = np.linspace(0, max_load, steps)

# Result containers
u2 = []
u3 = []
N_force = []
applied_loads = []

# Simulation loop
for f_applied in f_array:
    force_vector = np.array([0, f_applied, 0, 0])
    K_BCs = K_global[1:3, 1:3]
    F_BCs = force_vector[1:3]

    # Initial approximation
    u_trial = np.linalg.solve(K_BCs, F_BCs)
    gap = initial_gap - (u_trial[0] - u_trial[1])

    if gap < 0:  # Contact detected
        def residuals(u_vals, f_vals, K_sub, penalty_coeff, gap_zero):
            u_a, u_b = u_vals
            f1, f2 = f_vals
            k11, k12 = K_sub[0]
            k21, k22 = K_sub[1]
            gap_now = gap_zero - (u_a - u_b)
            res1 = f1 - k11 * u_a - k12 * u_b + penalty_coeff * gap_now
            res2 = f2 - k21 * u_a - k22 * u_b - penalty_coeff * gap_now
            return [res1, res2]

        u_cont = fsolve(residuals, u_trial, args=(F_BCs, K_BCs, penalty, initial_gap))
        final_gap = initial_gap - (u_cont[0] - u_cont[1])
        N = -penalty * final_gap
        u_final = u_cont
    else:
        N = 0
        u_final = u_trial

    # Store results
    u2.append(u_final[0])
    u3.append(-u_final[1])
    N_force.append(N)
    applied_loads.append(f_applied)

# Compute axial stresses
stress_1 = E * np.array(u2) / l1
stress_2 = E * np.array(u3) / l2

# Output final values
print("CONT FORCE: ", N_force[-1])
print("DISP B1: ", u2[-1])
print("DISP B2: ", u3[-1])
print("STRESS B1: ", stress_1[-1])
print("STRESS B2: ", stress_2[-1])

# Plot F - u
plt.figure(figsize=(8, 5))
plt.plot(u2, applied_loads)
plt.xlabel('Displacement at node 2 (mm)')
plt.ylabel('External Force P (N)')
plt.grid()
plt.show()

# Plot F - R
plt.figure(figsize=(8, 5))
plt.plot(N_force, applied_loads, color='orange')
plt.xlabel('Contact Force (N)')
plt.ylabel('Force P (N)')
plt.grid()
plt.show()

# Plot stresses vs force
plt.figure(figsize=(8, 5))
plt.plot(applied_loads, stress_1, label='Beam 1')
plt.plot(applied_loads, stress_2, label='Beam 2')
plt.xlabel('Contact Force (N)')
plt.ylabel('Principal Stress (MPa)')
plt.grid()
plt.legend()
plt.show()
