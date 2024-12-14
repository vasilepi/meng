import matplotlib.pyplot as plt
from bonus_implicit import bon_implicit
from bonus_explicit import bon_explicit
import numpy as np

# Get results
x_implicit, y_implicit, bl_implicit, d1_implicit, d2_implicit, cf_implicit, u_implicit, v_implicit = bon_implicit()
x_explicit, y_explicit, bl_explicit, d1_explicit, d2_explicit, cf_explicit, u_explicit, v_explicit = bon_explicit()

# grid
X_explicit, Y_explicit = np.meshgrid(x_explicit, y_explicit)
X_implicit, Y_implicit = np.meshgrid(x_implicit, y_implicit)



# Blasius
L = 8
Vinf = 1
nu = 1.5e-5
x_blasius = np.linspace(0.01, L, 16000)
Rex = Vinf * x_blasius / nu
bl_blasius = 5 * x_blasius / np.sqrt(Rex)
d1_blasius = 1.72 * x_blasius / np.sqrt(Vinf * x_blasius / nu + np.finfo(float).eps)
d2_blasius = 0.664 * x_blasius / np.sqrt(Vinf * x_blasius / nu + np.finfo(float).eps)
cf_blasius = 0.664 / np.sqrt(Vinf * x_blasius / nu + np.finfo(float).eps)
# ALL PLOTS
#grid vis
plt.figure()
plt.plot(X_explicit, Y_explicit, color='black')  # Vertical lines
plt.plot(X_explicit.T, Y_explicit.T, color='black')  # Horizontal lines
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.title('Visualization of the Physical Grid')
plt.grid(False)
plt.ylim(0, 0.1)
plt.xlim(0,0.03)
plt.show()

#delta
plt.figure()
plt.plot(x_explicit, bl_explicit, color='green', label="CFD Boundary Layer Thickness")
plt.plot(x_blasius, bl_blasius, color="blue", label="Blasius Boundary Layer Thickness")
plt.xlabel("x (m)")
plt.ylabel("Boundary Layer Thickness δ (m)")
plt.title("Boundary Layer Thickness over Flat Plate [Explicit]")
plt.legend()
plt.grid(True)
plt.show()

plt.figure()
plt.plot(x_implicit, bl_implicit, color='red', label="CFD Boundary Layer Thickness")
plt.plot(x_blasius, bl_blasius, color="blue", label="Blasius Boundary Layer Thickness")
plt.xlabel("x (m)")
plt.ylabel("Boundary Layer Thickness δ (m)")
plt.title("Boundary Layer Thickness over Flat Plate [Implicit]")
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 5))
plt.plot(x_implicit, bl_implicit, label="Implicit CFD", color="red")
plt.plot(x_explicit, bl_explicit, label="Explicit CFD", color="green")
plt.plot(x_blasius, bl_blasius, label="Blasius", color="blue")
plt.title("Boundary Layer Thickness over Flat Plate Comparison")
plt.xlabel("x (m)")
plt.ylabel("Boundary Layer Thickness δ (m)")
plt.legend()
plt.grid(True)
plt.show()

#delta1
plt.figure(figsize=(10, 5))
plt.plot(x_implicit, d1_implicit, label="Implicit CFD", color="red")
plt.plot(x_explicit, d1_explicit, label="Explicit CFD", color="green")
plt.plot(x_blasius, d1_blasius, label="Blasius", color="blue")
plt.title("δ1 Comparison")
plt.xlabel("x (m)")
plt.ylabel("Boundary Layer Thickness δ1 (m)")
plt.legend()
plt.grid(True)
# plt.ylim(0,0.02)
# plt.xlim(0,4)
plt.show()


#delta2
plt.figure(figsize=(10, 5))
plt.plot(x_implicit, d2_implicit, label="Implicit CFD", color="red")
plt.plot(x_explicit, d2_explicit, label="Explicit CFD", color="green")
plt.plot(x_blasius, d2_blasius, label="Blasius", color="blue")
plt.title("δ2 Comparison")
plt.xlabel("x (m)")
plt.ylabel("Boundary Layer Thickness δ2 (m)")
plt.legend()
plt.grid(True)
# plt.ylim(0,0.004)
# plt.xlim(0,2)
plt.show()

#cf
plt.figure(figsize=(10, 5))
plt.plot(x_implicit, cf_implicit, label="Implicit CFD", color="red")
plt.plot(x_explicit, cf_explicit, label="Explicit CFD", color="green")
plt.plot(x_blasius, cf_blasius, label="Blasius", color="blue")
plt.title("Friction Coefficient Comparison")
plt.xlabel("x (m)")
plt.ylabel("Friction Coefficient cf")
plt.legend()
plt.grid(True)
# plt.ylim(0,0.02)
# plt.xlim(0,2)
plt.show()
