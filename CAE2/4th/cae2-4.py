import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import StateSpace, lsim
from scipy.linalg import eigh, inv, solve

# Sampling
W1 = 10
W2 = 15
fs = 100 * max(W1, W2)
dt = 1 / fs
t = np.arange(0, 1 + dt, dt)

# Load
r1disp = 0.1 * np.sin(2*np.pi*W1 * t)
r2disp = 0.15 * np.sin(2*np.pi*W2 * t)
plt.figure()
plt.subplot(2, 1, 1)
plt.plot(t, r1disp, 'b', linewidth=1)
plt.ylabel('F (N)')
plt.xlabel('Time (sec)')
plt.title('DOF 1 Excitation')

plt.subplot(2, 1, 2)
plt.plot(t, r2disp, 'b', linewidth=1)
plt.ylabel('F (N)')
plt.xlabel('Time (sec)')
plt.title('DOF 2 Excitation')
plt.tight_layout()
plt.show()

# Model parameters
L = 4.5
m = 1500
Ig = 1805
m1 = 105
m2 = 110

k1 = 24.2E3
k2 = 27.1E3
kw1 = 360E3
kw2 = 340E3

c1 = 3E3
c2 = 2.8E3
cw1 = 0.15E3
cw2 = 0.16E3

ro1 = 0.100
ro2 = 0.150
R1 = cw1 * 2* np.pi*W1* ro1 * np.cos(2*np.pi*W1 * t) + kw1 * ro1 * np.sin(2*np.pi*W1 * t)
R2 = kw2 * ro2 * np.sin(2*np.pi*W2 * t) + cw2 * ro2 *2*np.pi*W2* np.cos(2*np.pi*W2 * t)
Rtot = np.vstack((R1, R2)).T

M = np.array([
    [m1, 0, 0, 0],
    [0, m2, 0, 0],
    [0, 0, m, 0],
    [0, 0, 0, Ig]
])

C = np.array([
    [c1 + cw1, 0, -c1, c1 * L / 3],
    [0, c2 + cw2, -c2, -2 * c2 * L / 3],
    [-c1, -c2, c1 + c2, (2 * c2 - c1) * L / 3],
    [c1 * L / 3, -2 * c2 * L / 3, (2 * c2 - c1) * L / 3, (4 * c2 + c1) * (L / 3) ** 2]
])

K = np.array([
    [k1 + kw1, 0, -k1, k1 * L / 3],
    [0, k2 + kw2, -k2, -2 * k2 * L / 3],
    [-k1, -k2, k1 + k2, (2 * k2 - k1) * L / 3],
    [k1 * L / 3, -2 * k2 * L / 3, (2 * k2 - k1) * L / 3, (4 * k2 + k1) * (L / 3) ** 2]
])

ndof = 4
next = 2

# Input location
lforce = [0, 1]  # zero-based indexing for Python
LN = np.zeros((ndof, next))
for i in range(ndof):
    for j in range(next):
        if lforce[j] == i:
            LN[i, j] = 1

Minv = inv(M)

AN = np.block([
    [-Minv @ C, -Minv @ K],
    [np.eye(ndof), np.zeros((ndof, ndof))]
])
BN = np.vstack([
    -Minv @ LN,
    np.zeros((ndof, next))
])
CN = np.block([
    [-Minv @ C, -Minv @ K],
    [np.eye(ndof), np.zeros((ndof, ndof))],
    [np.zeros((ndof, ndof)), np.eye(ndof)]
])
DN = np.vstack([
    Minv @ LN,
    np.zeros((ndof, next)),
    np.zeros((ndof, next))
])

# State-space system
half_car = StateSpace(AN, BN, CN, DN)

# Modal analysis
eigvals, Phi = eigh(K, M)
freq = np.sqrt(np.real(eigvals))
indx = np.argsort(freq)
modfreq = freq[indx] / (2 * np.pi)

PhiSorted = Phi[:, indx]
num = PhiSorted.T @ C @ PhiSorted
den = PhiSorted.T @ K @ PhiSorted
zeta = np.diag(1/2/freq[indx] * num/den)

print(zeta)

freq_sorted = freq[indx]
omega1 = freq_sorted[0]
omega2 = freq_sorted[1]
z1 = zeta[0]
z2 = zeta[1]
A = np.array([[1/(2*omega1), omega1/2],
              [1/(2*omega2), omega2/2]])
b = np.array([z1, z2])
alpha, beta = np.linalg.solve(A, b)


print("Rayleig a:\n", alpha)
print("Rayleig b:\n", beta)



print("Normal Modes (eigenvectors):\n", Phi[:, indx])
print("Modal Frequencies (Hz):\n", modfreq)
print("Modal Frequencies (rad/s):\n", freq)

# Simulation
SPC = np.zeros((2 * ndof,))
t_out, y_out, x_out = lsim(half_car, Rtot, t, X0=SPC)
y_out[:, 0] *= -1  # Front Wheel accel
y_out[:, 1] *= -1  # Rear Wheel accel
y_out[:, 2] *= -1  # Pitch angular accel

# Plotting simulation results
labels = [
    'Wheel Measurement DOF 1',
    'Wheel Measurement DOF 2',
    'Translational Body DOF 3',
    'Rotational Body DOF 4'
]
colors = ['r', 'g', 'k', 'm']

plt.figure(figsize=(10, 8))
for i in range(4):
    plt.subplot(2, 2, i+1)
    plt.plot(t_out, y_out[:, i], colors[i], linewidth=1)
    plt.title(labels[i], fontsize=12, fontweight='bold', fontname='Times New Roman')
    plt.ylabel('Acceleration [m/s²]', fontsize=12, fontweight='bold', fontname='Times New Roman')
    plt.xlabel('Time [s]')
    plt.grid(True)
    plt.legend(['Linear-SS'])
    plt.tick_params(axis='both', which='major', labelsize=10)
plt.tight_layout()
plt.show()

displacements = x_out[:, 4:]  # Columns x1 to x4

# Apply sign correction
displacements[:, 2] *= -1
displacements[:, 1] *= -1
displacements[:, 0] *= -1
disp_titles = [
    "Front Wheel Displacement (DOF 1)",
    "Rear Wheel Displacement (DOF 2)",
    "Body Vertical Displacement (DOF 3)",
    "Body Pitch Angle θ (DOF 4)"
]
disp_ylabels = [
    "Displacement [m]",
    "Displacement [m]",
    "Displacement [m]",
    "Pitch Angle [rad]"
]

plt.figure(figsize=(10, 8))
for i in range(4):
    plt.subplot(2, 2, i+1)
    plt.plot(t_out, displacements[:, i], linewidth=1)
    plt.title(disp_titles[i], fontsize=12, fontweight='bold')
    plt.ylabel(disp_ylabels[i], fontsize=12)
    plt.xlabel("Time [s]", fontsize=10)
    plt.grid(True)
plt.tight_layout()
plt.show()

np.set_printoptions(threshold=np.inf)  # Show full array without summarizing

# print(t_out)
# print(y_out[:, 3])