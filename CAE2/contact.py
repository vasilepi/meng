# SCRIPT FOR 2nd assignement on CAE2

##### PART B THEORETICAL
### HERTZ LINEAR

import numpy as np
from matplotlib import pyplot as plt

r1 = 5 # mm
r2 = 5 # mm
l = 10 # mm
f = 1000 # N

n1 = 0.3
n2 = 0.3
E1 = 210000 # MPa
E2 = 210000 # MPa


E_star = 1/((1-n1**2)/E1 + (1-n2**2)/E2)
R_star = 1/(1/r1 + 1/r2)

# print(E_star, R_star)
alin = np.sqrt(4*f*R_star/np.pi/E_star/l)
# print(alin)

polin = 2*f/np.pi/alin/l
xlin = np.linspace(-alin, alin)
plin = polin * np.sqrt(1-(xlin/alin)**2)

ulin = np.linspace(0, 10)
Flin = np.pi/4 * E_star*l * ulin

zlin = np.linspace(-10,10)
slinx1 = -2* n1 * polin * (np.sqrt(1+(zlin/alin)**2)- np.abs(zlin/alin))
slinx2 = -2* n2 * polin * (np.sqrt(1+(zlin/alin)**2)- np.abs(zlin/alin))
sliny = -polin * ((1+2*(zlin/alin)**2)/np.sqrt(1+(zlin/alin)**2) - 2* np.abs(zlin/alin))
slinz = -polin/np.sqrt(1+(zlin/alin)**2)



# plot
plt.figure(figsize=(9,5))
plt.grid(True, which='both')
plt.plot(xlin, plin)
plt.legend()
plt.xlabel(r"$x\; (mm)$")
plt.ylabel(r"$P\; (MPa)$")
plt.show()

plt.figure(figsize=(9,5))
plt.grid(True, which='both')
plt.plot(ulin, Flin)
plt.legend()
plt.xlabel(r"$u\; (mm)$")
plt.ylabel(r"$F\; (N)$")
plt.show()

plt.figure(figsize=(9,5))
plt.grid(True, which='both')
plt.plot(zlin, slinx1)
plt.legend()
plt.xlabel(r"$z\; (mm)$")
plt.ylabel(r"$sx\; (MPa)$")
plt.show()

plt.figure(figsize=(9,5))
plt.grid(True, which='both')
plt.plot(zlin, sliny)
plt.legend()
plt.xlabel(r"$z\; (mm)$")
plt.ylabel(r"$sy\; (MPa)$")
plt.show()

plt.figure(figsize=(9,5))
plt.grid(True, which='both')
plt.plot(zlin, slinz)
plt.legend()
plt.xlabel(r"$z\; (mm)$")
plt.ylabel(r"$sz\; (MPa)$")
plt.show()
# point

r1p = 5 # mm
fp = 1000 # N

E_starp = 1/((1-n1**2)/E1 + (1-n2**2)/E2)
R_starp = 1/(1/r1p)

# print(E_star, R_star)
ap = (3*fp*R_starp/E_starp/4)
print(ap)

pop = 3*fp/2/np.pi/ap**2
rp = np.linspace(-ap, ap)
pp = pop * np.sqrt(1-(rp/ap)**2)

up = np.linspace(0, 10)
Fp = 4*np.pi/3 * E_starp*R_starp**(0.5) * up**(3/2)

zp = np.linspace(-10,10)
sp1 = -pop * ((1-np.abs(zp/ap)*np.atan(1/np.abs(zp/ap)))*(1+n1)-1/(2*(1+(zp/ap)**2)))
sp2 = -pop * ((1-np.abs(zp/ap)*np.atan(1/np.abs(zp/ap)))*(1+n2)-1/(2*(1+(zp/ap)**2)))
spz = -pop/(1+(zp/ap)**2)



# plot
# plt.figure(figsize=(9,5))
# plt.grid(True, which='both')
# plt.plot(rp, pp)
# plt.legend()
# plt.xlabel(r"$x\; (mm)$")
# plt.ylabel(r"$P\; (MPa)$")
# plt.show()

# plt.figure(figsize=(9,5))
# plt.grid(True, which='both')
# plt.plot(up, Fp)
# plt.legend()
# plt.xlabel(r"$u\; (mm)$")
# plt.ylabel(r"$F\; (N)$")
# plt.show()

# plt.figure(figsize=(9,5))
# plt.grid(True, which='both')
# plt.plot(zp, sp1)
# plt.legend()
# plt.xlabel(r"$z\; (mm)$")
# plt.ylabel(r"$sx\; (MPa)$")
# plt.show()

# plt.figure(figsize=(9,5))
# plt.grid(True, which='both')
# plt.plot(zp, spz)
# plt.legend()
# plt.xlabel(r"$z\; (mm)$")
# plt.ylabel(r"$sz\; (MPa)$")
# plt.show()