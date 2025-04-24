# weld 1st ass
import numpy as np
from matplotlib import pyplot as plt

#kN
Fm = np.array([28.2, 37.6, 34.2, 31.2, 40.2, 22.1, 18.1, 20.1, 21.1, 26.2, 16.5, 23.2, 24.2])
Fa = np.array([27.9, 37.3, 33.9, 30.9, 39.9, 22, 18, 20, 21, 26, 16.4, 23, 24])


N = np.array([89000, 29000, 41000, 73000, 23000, 190000, 603000, 453000, 575000, 167000, 1470000, 335000, 246000])


### NOMINAL

Fmin = Fm - Fa
Fmax = Fm + Fa

t = 2.24 # mm
l = 50 # mm
aw = 1.6 # mm
A = 223.983 # mm^2

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



### STRUCTURAL HOT SPOT 

## measured at 0.4t and 1.0t -> shs = [s_at_0.4t, s_at_1.0t]
## shs = 1.67*s_0.4t - 0.67*s_1.0t

# normal x

shs100_comp = np.array([9.679E-2, 2.892E-2])    # 0
shs200_comp = np.array([1.936E-1, 5.785E-2])    # 1
shs300_comp = np.array([2.904E-1, 8.677E-2])    # 2
shs32900_comp = np.array([3.184e+1, 9.516])     # 3
shs36100_comp = np.array([3.494E+1, 1.044E+1])  # 4
shs40100_comp = np.array([3.881E+1, 1.16E+1])   # 5
shs42100_comp = np.array([4.075E+1, 1.218E+1])  # 6
shs44100_comp = np.array([4.268E+1, 1.276E+1])  # 7
shs46200_comp = np.array([4.472E+1, 1.336E+1])  # 8
shs48200_comp = np.array([4.665E+1, 1.394E+1])  # 9
shs52200_comp = np.array([5.052E+1, 1.510E+1])  # 10
shs56100_comp = np.array([5.43E+1, 1.623E+1])   # 11
shs62100_comp = np.array([6.01E+1, 1.796E+1])   # 12
shs68100_comp = np.array([6.591E+1, 1.97E+1])   # 13
shs74900_comp = np.array([7.249E+1, 2.166E+1])  # 14
shs80100_comp = np.array([7.753E+1, 2.317E+1])  # 15

shs_full = np.array((shs100_comp,shs200_comp,shs300_comp,shs32900_comp,shs36100_comp,shs40100_comp,shs42100_comp,shs44100_comp,shs46200_comp,shs48200_comp,shs52200_comp,shs56100_comp,shs62100_comp,shs68100_comp,shs74900_comp,shs80100_comp))

shs=[]
for i in range(0,len(shs_full)):
    shs_extrap = 1.67*shs_full[i,0]-0.67*shs_full[i,1]
    shs.append(shs_extrap)

shs_idx = np.array(([2,11],[2,14],[2,13],[2,12],[2,15],[0,7],[0,4],[0,5],[0,6],[1,10],[0,3],[1,8],[1,9]))

Ds_shs = []
for j in range(0,len(shs_idx)):
    Ds_shs.append(shs[shs_idx[j,1]]-shs[shs_idx[j,0]]) 



# major principal
shs100_comp2 = np.array([7.342e-1, 6.038e-1])    # 0
shs200_comp2 = 2*shs100_comp2 # 1
shs300_comp2 = 3*shs100_comp2    # 2
shs32900_comp2 = 329*shs100_comp2     # 3
shs36100_comp2 = 361*shs100_comp2  # 4
shs40100_comp2 = 401*shs100_comp2   # 5
shs42100_comp2 = 421*shs100_comp2  # 6
shs44100_comp2 = 441*shs100_comp2  # 7
shs46200_comp2 = 462*shs100_comp2  # 8
shs48200_comp2 = 482*shs100_comp2  # 9
shs52200_comp2 = 522*shs100_comp2  # 10
shs56100_comp2 = 561*shs100_comp2   # 11
shs62100_comp2 = 621*shs100_comp2   # 12
shs68100_comp2 = 681*shs100_comp2   # 13
shs74900_comp2 = 749*shs100_comp2  # 14
shs80100_comp2 = 801*shs100_comp2  # 15

shs_full2 = np.array((shs100_comp2,shs200_comp2,shs300_comp2,shs32900_comp2,shs36100_comp2,shs40100_comp2,shs42100_comp2,shs44100_comp2,shs46200_comp2,shs48200_comp2,shs52200_comp2,shs56100_comp2,shs62100_comp2,shs68100_comp2,shs74900_comp2,shs80100_comp2))

shs2=[]
for i in range(0,len(shs_full2)):
    shs_extrap2 = 1.67*shs_full2[i,0]-0.67*shs_full2[i,1]
    shs2.append(shs_extrap2)

shs_idx2 = np.array(([2,11],[2,14],[2,13],[2,12],[2,15],[0,7],[0,4],[0,5],[0,6],[1,10],[0,3],[1,8],[1,9]))

Ds_shs2 = []
for j in range(0,len(shs_idx2)):
    Ds_shs2.append(shs2[shs_idx2[j,1]]-shs2[shs_idx2[j,0]]) 


# # statistics
x_hat_shs = np.log10(Ds_shs2)
nn_shs = len(x_hat_shs)
#statistical
x_shs_ = np.mean(x_hat_shs)

x2_shs = x_hat_shs**2

xy_shs = x_hat_shs*y_hat

# y^
b_shs = (sum(xy_shs)-sum(x_hat_shs)*sum(y_hat)/nn_shs)/(sum(x2_shs)-nn_shs*x_shs_**2)
a_shs = sum(y_hat)/nn_shs - b_shs * sum(x_hat_shs)/nn_shs

y_pred_shs = a_shs + b_shs*x_hat_shs
N_pred_shs = 10**y_pred_shs

yx_shs_ = y_hat + b_shs*(x_shs_ - x_hat_shs)
yx_y2_shs = (yx_shs_ - y_)**2

sigma_shs = np.sqrt(sum(yx_y2_shs)/nn_shs)

# 97.7%
y977_shs = y_pred_shs - 2.27*sigma_shs
N977_shs = 10**(y977_shs)






### FAT


## Stress-Ratio factor
R = 0 # give


if R <-1:
    fR = 1.6
elif R>0.5:
    fR = 1
else:
    fR = -0.4*R + 1.2

## Wall thickness factor
# give n
n_thick = 0.3

if l/t<2:
    teff=t
else:
    teff=np.max([0.5*l,t])

# print(teff)

tref = 25
ft = (tref/teff)**n_thick


# print(ft*fR)

#FAT71 nominal
Ds_fat = 71 # MPa
m = 3 # Steel

Cfat = 2E6 * Ds_fat**m
fat71log = np.linspace(1e4,2e6)
fat71 = (Cfat/fat71log)**(1/m) 
fat71mod = (Cfat/fat71log)**(1/m) * ft *fR

#FATx hot-spot

# ref-assess
# at 100N
# shs_ref = 0.057
# shs_assess = shs2[0]

# # print(shs2[0])

# Ds_shs_fat = shs_ref/shs_assess * 100 # MPa
# Cfat_shs = 2E6 * Ds_shs_fat**m
# fatHSlog = np.linspace(1e4,2e6)
# fatHS = (Cfat_shs/fatHSlog)**(1/m) * ft * fR

C100 = 2E6 * 100**m
fat100log = np.linspace(1e4,2e6)
fat100mod = (C100/fat100log)**(1/m) * ft * fR
fat100 = (C100/fat100log)**(1/m)


# PLOT
# plt.figure(figsize=(8,5))
# plt.scatter(y_hat,x_hat)
# plt.plot(y_pred, x_hat)
# plt.plot(y977,x_hat)
# plt.grid()
# plt.show()

plt.figure(figsize=(9,5))
plt.grid(True, which='both')
plt.scatter(N, Ds, color = 'black')
plt.plot(N_pred,Ds)
plt.plot(N977, Ds)
plt.plot(fat71log,fat71)
plt.plot(fat71log,fat71mod)
plt.xscale('log')
plt.yscale('log')
plt.legend(['Experimental Data', '50% Probability', '97.7% Probability', 'FAT71', 'FAT71 - Modified'])
plt.xlabel(r"$N$")
plt.ylabel(r"$\Delta \sigma$ (MPa)")
plt.show()

plt.figure(figsize=(9,5))
plt.grid(True, which='both')
# plt.scatter(N, Ds_shs, color='black', label = "Normal to weld toe")
plt.scatter(N, Ds_shs2, color='black', label = "Major Principal")
plt.plot(N_pred_shs,Ds_shs2, label = "50% Probability")
plt.plot(N977_shs, Ds_shs2, label = "97.7% Probability")
plt.plot(fat100log,fat100, label = "FAT100")
plt.plot(fat100log,fat100mod, label = "FAT100 - Modified")
# plt.plot(fatHSlog,fatHS, label = "FAT" + str(round(Ds_shs_fat,2)))
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r"$N$")
plt.ylabel(r"$\Delta \sigma$ (MPa)")
plt.show()