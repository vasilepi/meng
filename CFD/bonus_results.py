import matplotlib.pyplot as plt
from bonus_implicit import bon_implicit
from bonus_explicit import bon_explicit
from flat_plate_bl_implicit import run_implicit
import numpy as np
from matplotlib import patheffects

# Get results
x_implicit, y_implicit, bl_implicit, d1_implicit, d2_implicit, cf_implicit, u_implicit, v_implicit = bon_implicit()
x_explicit, y_explicit, bl_explicit, d1_explicit, d2_explicit, cf_explicit, u_explicit, v_explicit = bon_explicit()
x1_implicit, y1_implicit, bl1_implicit, d11_implicit, d21_implicit, cf1_implicit, u1_implicit, v1_implicit = run_implicit()

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
plt.ylim(0, 0.08)
plt.xlim(0,0.03)
plt.show()

#delta
# plt.figure()
# plt.plot(x_implicit, bl_implicit, color='green', label=r"Non-constant $\Delta y$")
# plt.plot(x1_implicit, bl1_implicit, color="blue", label="Constant $\Delta y$")
# plt.xlabel("x (m)")
# plt.ylabel("Boundary Layer Thickness δ (m)")
# plt.title("Boundary Layer Thickness over Flat Plate Comparison")
# plt.legend()
# plt.grid(True)
# plt.show()

plt.figure()
plt.plot(x_implicit, bl_implicit, color='red', label="CFD Boundary Layer Thickness")
plt.plot(x_blasius, bl_blasius, color="blue", label="Blasius Boundary Layer Thickness")
plt.xlabel("x (m)")
plt.ylabel("Boundary Layer Thickness δ (m)")
plt.title("Boundary Layer Thickness over Flat Plate")
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 5))
plt.plot(x_implicit, bl_implicit, color='red', label=r"Non-constant $\Delta y$")
plt.plot(x1_implicit, bl1_implicit, color="green", label=r"Constant $\Delta y$")
plt.plot(x_blasius, bl_blasius, label="Blasius", color="blue")
plt.title("Boundary Layer Thickness over Flat Plate Comparison")
plt.xlabel("x (m)")
plt.ylabel("Boundary Layer Thickness δ (m)")
plt.legend()
plt.grid(True)
plt.show()

#delta1
plt.figure(figsize=(10, 5))
plt.plot(x_implicit, d1_implicit, label=r"Non-constant $\Delta y$", color="red")
plt.plot(x1_implicit, d11_implicit, label=r"Constant $\Delta y$", color="green")
plt.plot(x_blasius, d1_blasius, label="Blasius", color="blue")
plt.title("δ1 Comparison")
plt.xlabel("x (m)")
plt.ylabel("Boundary Layer Thickness δ1 (m)")
plt.legend()
plt.grid(True)
# plt.ylim(0,0.02)
# plt.xlim(0,4)
plt.show()

# #delta1     ZOOM
# plt.figure(figsize=(10, 5))
# plt.plot(x_implicit, d1_implicit, label=r"Non-constant $\Delta y$", color="red")
# plt.plot(x1_implicit, d11_implicit, label=r"Constant $\Delta y$", color="green")
# plt.plot(x_blasius, d1_blasius, label="Blasius", color="blue")
# plt.title("δ1 Comparison")
# plt.xlabel("x (m)")
# plt.ylabel("Boundary Layer Thickness δ1 (m)")
# plt.legend()
# plt.grid(True)
# plt.ylim(0,0.02)
# plt.xlim(0,4)
# plt.show()


#delta2
plt.figure(figsize=(10, 5))
plt.plot(x_implicit, d2_implicit, label=r"Non-constant $\Delta y$", color="red")
plt.plot(x1_implicit, d21_implicit, label=r"Constant $\Delta y$", color="green")
plt.plot(x_blasius, d2_blasius, label="Blasius", color="blue")
plt.title("δ2 Comparison")
plt.xlabel("x (m)")
plt.ylabel("Boundary Layer Thickness δ2 (m)")
plt.legend()
plt.grid(True)
# plt.ylim(0,0.004)
# plt.xlim(0,2)
plt.show()

#delta2 ZOOM
red_effect = [patheffects.withTickedStroke(spacing=30, angle=135, length=0.2)]
green_effect = [patheffects.withTickedStroke(spacing=20, angle=215, length=0.4)]

plt.figure(figsize=(10, 5))
line11, = plt.plot(x_implicit, d2_implicit, label=r"Non-constant $\Delta y$", path_effects=red_effect, color="red")
line12, = plt.plot(x1_implicit, d21_implicit, label=r"Constant $\Delta y$", path_effects=green_effect, color="green")
plt.plot(x_blasius, d2_blasius, label="Blasius", color="blue")
plt.title("δ2 Comparison")
plt.xlabel("x (m)")
plt.ylabel("Boundary Layer Thickness δ2 (m)")
legend1 = plt.legend()
for text, line in zip(legend1.get_texts(), [line11, line12]):
    text.set_path_effects(line.get_path_effects())
plt.grid(True)
plt.ylim(0, 0.004)
plt.xlim(0, 2)
plt.show()

#cf
plt.figure(figsize=(10, 5))
plt.plot(x_implicit, cf_implicit, label=r"Non-constant $\Delta y$", color="red")
plt.plot(x1_implicit, cf1_implicit, label=r"Constant $\Delta y$", color="green")
plt.plot(x_blasius, cf_blasius, label="Blasius", color="blue")
plt.title("Friction Coefficient Comparison")
plt.xlabel("x (m)")
plt.ylabel("Friction Coefficient cf")
plt.legend()
plt.grid(True)
# plt.ylim(0,0.02)
# plt.xlim(0,2)
plt.show()

#cf         ZOOM
plt.figure(figsize=(10, 5))
line12, = plt.plot(x_implicit, cf_implicit, label=r"Non-constant $\Delta y$", path_effects=red_effect, color="red")
line22, = plt.plot(x1_implicit, cf1_implicit, label=r"Constant $\Delta y$", path_effects=green_effect, color="green")
plt.plot(x_blasius, cf_blasius, label="Blasius", color="blue")
plt.title("Friction Coefficient Comparison")
plt.xlabel("x (m)")
plt.ylabel("Friction Coefficient cf")
legend2 = plt.legend()
for text, line in zip(legend2.get_texts(), [line12, line22]):
    text.set_path_effects(line.get_path_effects())
plt.grid(True)
plt.ylim(0,0.01)
plt.xlim(0.2,1.5)
plt.show()

