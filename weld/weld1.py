# weld 1st ass
import numpy as np
from matplotlib import pyplot as plt

#kN
Fm = np.array([28.2, 37.6, 34.2, 31.2, 40.2, 22.1, 18.1, 20.1, 21.1, 26.2, 16.5, 23.2, 24.2])
Fa = np.array([27.9, 37.3, 33.9, 30.9, 39.9, 22, 18, 20, 21, 26, 16.4, 23, 24])


N = np.array([89000, 29000, 41000, 73000, 23000, 190000, 603000, 453000, 575000, 167000, 1470000, 335000, 246000])


#Nominal

Fmin = Fm - Fa
Fmax = Fm + Fa


l = 50 # mm
aw = 1.6 # mm
A = l*aw # mm^2

smin = Fmin/A * 1000 # MPa
smax = Fmax/A *1000 # MPa

Ds = smax - smin


x_hat = np.log10(Ds)
y_hat = np.log10(N)
nn = len(x_hat)
#statistical
x_ = np.mean(x_hat)
y_ = np.mean(y_hat)

x2 = x_hat**2
y2 = y_hat**2
xy = x_hat*y_hat

# y^
b = (sum(xy)-sum(x_hat)*sum(y_hat)/nn)/(sum(x2)-nn*x_**2)
a = sum(y_hat)/nn -b * sum(x_hat)/nn

y_pred = a + b*x_hat
N_pred = 10**y_pred

yx_ = y_hat + b*(x_ - x_hat)
yx_y2 = (yx_ - y_)**2

sigma = np.sqrt(sum(yx_y2)/nn)

# 97.7%
y977 = y_pred - 2.27*sigma
N977 = 10**(y977)


# FAT

#FAT71 nominal
Ds_fat = 71 # MPa
m = 3 # Steel

Cfat = 2E6 * Ds_fat**m
fatlog = np.linspace(1e4,1e7)
fat71 = (Cfat/fatlog)**(1/m)



# PLOT
# plt.figure(figsize=(8,5))
# plt.scatter(y_hat,x_hat)
# plt.plot(y_pred, x_hat)
# plt.plot(y977,x_hat)
# plt.grid()
# plt.show()

plt.figure(figsize=(8,5))
plt.grid(True, which='both')
plt.scatter(N, Ds, color = 'black')
plt.plot(N_pred,Ds)
plt.plot(N977, Ds)
plt.plot(fatlog,fat71)
plt.xscale('log')
plt.yscale('log')
plt.legend(['Experimental Data', '50% Probability', '97.7% Probability', 'FAT71'])
plt.xlabel(r"$N$")
plt.ylabel(r"$\Delta \sigma$ (MPa)")
plt.show()