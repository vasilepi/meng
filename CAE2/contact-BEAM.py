from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import minimize

# Constants and properties
Farr = np.linspace(0, 6000, 100)  # Applied forces [N]
A1 = 1.8  # Area bar 1 [mm^2]
A2 = 3.1  # Area bar 2 [mm^2]
l1 = 500  # Length bar 1 [mm]
l2 = 1000 # Length bar 2 [mm]
l = 1505  # Distance between walls [mm]
E = 210e3 # Young's modulus [MPa]
e = 1e8   # Penalty stiffness [N/mm]

# Stiffnesses
k1 = E * A1 / l1
k2 = E * A2 / l2

# Lists to store results
u_vals = []
uu_vals = []
s1_vals = []
s2_vals = []
N_vals = []

# Loop over applied force
for f in Farr:
    # Free displacement of bar 1
    u = f / k1
    gap = l - (l1 + l2 + u)

    if gap > 0:
        # No contact
        N = 0
        s1 = f / A1
        s2 = 0
        u1 = u
        u2 = 0
    else:
        # Contact occurs
        u1 = minimize(lambda u: 0.5*k1*u**2 + 0.5*k2*(u - (l-l1-l2))**2 - f*u + 0.5*e*gap**2, u1).x[0] 
        u2 = (l1+l2+u1) - l
        gap = l - (l1+l2+u1)

        N = - e * gap
        s1 = (N_vals[-2]-f) / A1
        s2 = -N_vals[-2] / A2
    print(gap)
    u_vals.append(u1)
    uu_vals.append(u2)
    s1_vals.append(s1)
    s2_vals.append(s2)
    N_vals.append(N)



print(s1_vals)
print("Contact force:",N_vals[-1])
print("Beam 1 deformation:",u_vals[-1])
print("Beam 1 stress:",s1_vals[-1])
print("Beam 2 deformation:",uu_vals[-1])
print("Beam 2 stress:",s2_vals[-1])


plt.figure(1)
plt.title("Beam 1 deformation")
plt.plot(u_vals,Farr, label = "Beam 1 Deformation")
plt.ylabel("Force (N)")
plt.xlabel("Deformation (mm)")
plt.show()

plt.figure(2)
plt.title("Beam Stresses")
plt.plot(Farr,s1_vals,label = "Beam 1 Stress")
plt.plot(Farr,s2_vals,label = "Beam 2 Stress")
plt.xlabel("Force (N)")
plt.ylabel("Stress (MPa)")
plt.legend(loc = "lower center")
plt.show()

plt.figure(3)
plt.title("Contact force")
plt.plot(N_vals,Farr, label = "Beam 1 Deformation")
plt.ylabel("Force (N)")
plt.xlabel("Contact force (N)")
plt.show()
