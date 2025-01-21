import strain_life_utils_revisti as sl
import numpy as np
from matplotlib import pyplot as plt
import time
from scipy.optimize import fsolve
start_time = time.time()



def solveNf(strain_amplitude,sf,ef,b,c,E):

    eq = lambda Nf: strain_amplitude - sf/E * (2*Nf)**b - ef * (2*Nf)**c

    Nf_init = 1e4

    Nf_final = fsolve(eq, Nf_init)

    return Nf_final[0]

def solve_morr(strain_amplitude,mean_stress,sf,ef,b,c,E):

    eq = lambda Nf: strain_amplitude - (sf-mean_stress)/E * (2*Nf)**b - ef * (2*Nf)**c

    Nf_init = 1e4

    Nf_final = fsolve(eq, Nf_init)

    return Nf_final[0]

def solve_pswt(strain_amplitude,max_stress,sf,ef,b,c,E):

    eq = lambda Nf: max_stress*strain_amplitude*E - sf**2 * (2*Nf)**(2*b) - sf * ef * E * (2*Nf)**(b+c)

    Nf_init = 1e4

    Nf_final = fsolve(eq, Nf_init)

    return Nf_final[0]





         #90% #50% #10%
E = 210000
Ku = 834 #922 #834 #775
nu = 0.15
Kt = 1.574
s_normalization = 174 #MPa

sf = 853  #774 #853 #940
ef = 0.707 #0.502 #0.707 #0.996
b = -0.0971
c = -0.619
scale = 1


df = np.loadtxt('Strain Gauge Console.txt')
sequence = np.array(df[1:]) * s_normalization/100 * scale * Kt

plt.figure(figsize=(10,6))
DS = sequence[0]
si = sl.solve_Neuber(DS, 0, 0, E, Ku, nu, 'R-O')
ei = sl.ramberg_osgood(si, E, Ku, nu)

sl.ramberg_osgood_plot(si, E, Ku, nu)

# initialize vectors
DeltaS = []
start = []
direction = []
Loops = []


for i, S in enumerate(sequence[1:]):

    r_o_int = 0
    closing = 0
    DS_temp = S - sequence[i]
    DS = abs(DS_temp)
    dir = DS/DS_temp
    DS_ = DS

    while len(DeltaS)!=0 and DS_ > DeltaS[-1]:
        closing = 1

        Loops.append ([DeltaS[-1], start[-1][0], start[-1][1], direction[-1]])
        Ds = sl.solve_Neuber(DeltaS[-1], ei, si, E, Ku, nu, 'Masing')
        De = sl.masing(Ds, E, Ku, nu)
        Loops.append (np.array([DeltaS[-1], start[-1][0]+Ds*direction[-1], start[-1][1]+De*direction[-1], dir]))
        if len(DeltaS) == 1:
            DS = DS_ - DeltaS[-1]
            DS_ = DeltaS[-1]

            r_o_int = 1
            si = start[-1][0]
            ei = start[-1][1]
            del DeltaS[-1], start[-1], direction[-1]
            break

        DS_ = DeltaS[-2] + (DS_ - DeltaS[-1])
        Ds = sl.solve_Neuber(DS_, ei, si, E, Ku, nu, 'Masing')
        De = sl.masing(Ds, E, Ku, nu)
        si = start[-2][0]
        ei = start[-2][1]

        del DeltaS[-2], start[-2], direction[-2]
        del DeltaS[-1], start[-1], direction[-1]


    if closing == 1:
        if r_o_int != 1:
            DeltaS.append(DS_)
            start.append(np.array([si, ei]))
            direction.append(dir)
            si = start[-1][0] + Ds * dir
            ei = start[-1][1] + De * dir

    if closing != 1:
        Ds = sl.solve_Neuber(DS, ei, si, E, Ku, nu, 'Masing')
        De = sl.masing(Ds, E, Ku, nu)
        DeltaS.append(DS)
        start.append(np.array([si, ei]))
        direction.append(dir)

        si = si + Ds * dir
        ei = ei + De * dir

    if r_o_int == 1:
        si = sl.solve_Neuber(DS, ei, si, E, Ku, nu, 'R-O')
        ei = sl.ramberg_osgood(si, E, Ku, nu)
        sl.ramberg_osgood_plot(si, E, Ku, nu)

        DeltaS = []
        start = []
        direction = []


damage, damage_morr, damage_pswt = 0, 0, 0
for [DS, stress_from, strain_from, dir] in Loops:
    Ds_closed = sl.solve_Neuber(DS, 0, 0, E, Ku, nu, 'Masing')
    De_closed = sl.masing(Ds_closed, E, Ku, nu)
    ea = De_closed/2
    sm = stress_from + Ds_closed/2 * dir
    smax = max(stress_from, stress_from+Ds_closed*dir)

    # No mean stress effect
    Nf = solveNf(ea,sf,ef,b,c,E)
    damage = damage + 1/Nf
    # Mean stress effect - solve_morr
    Nf_morr = solve_morr(ea,sm,sf,ef,b,c,E)
    damage_morr = damage_morr + 1/Nf_morr
    # Mean stress effect - solve_pswt
    if sm > 0:
        Nf_pswt = solve_pswt(ea,smax,sf,ef,b,c,E)
        damage_pswt = damage_pswt + 1/Nf_pswt

CumDam = 1/damage
CumDam_Morr = 1/damage_morr
CumDam_pswt = 1/damage_pswt

print(CumDam, CumDam_Morr, CumDam_pswt)

Loops_cum = [[DeltaS[i], start[i][0], start[i][1], direction[i]] for i in range(len(DeltaS))] + Loops
for [DS, stress_from, strain_from, dir] in Loops_cum:
    Ds = sl.solve_Neuber(DS, strain_from, stress_from, E, Ku, nu, 'Masing')
    s = stress_from + Ds * dir
    sl.masing_plot(stress_from, strain_from, s, E, Ku, nu, dir)


plt.grid(True)
plt.title('Hysterisis Loops - Stress Strain')
plt.ylabel('Stress')
plt.xlabel('Strain')
plt.show()


# strain life
Nplot = np.arange(1e3, 1e8, 100)
eaplot = (sf/E * (2*Nplot)**b + ef*(2*Nplot)**c)*100 # %

# strain-life graph
plt.figure(figsize = (10,6))
plt.plot(Nplot, eaplot)
plt.title('Strain Life - 50%')
plt.ylabel('Strain %')
plt.xlabel('Cycles to failure')
plt.grid()
plt.xscale('log')
plt.ylim((0, 1))
plt.show()



# plot sequence
# t = np.linspace(sequence[0], np.max(sequence), len(sequence))
# plt.figure(figsize=(10,6))
# plt.plot(t, sequence, label = 'Elastic notch stress')
# plt.plot(t, sequence/Kt, label = 'Elastic stress')
# plt.title('Stress history')
# plt.ylabel('Stress %')
# plt.xlabel('Points (not physical time)')
# plt.legend()
# plt.grid(True)
# plt.show()


end_time = time.time()
elapsed_time = end_time - start_time
print(f"Elapsed time: {elapsed_time:.2f} seconds")

