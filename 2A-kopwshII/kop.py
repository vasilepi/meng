import strain_life_utils as sl
import numpy as np
from matplotlib import pyplot as plt
import time
import scipy.optimize as opt
from scipy.optimize import fsolve

start_time = time.time()

E = 210000
Ku = 922
nu = 0.15
b = -0.0971
c = -0.619
sf = 774
ef = 0.502
Kt = 1.574
s_normalization = 174
factor = 1
sequence = np.loadtxt('Strain Gauge Console.txt')
sequence = (np.array(sequence[:100])) * (s_normalization / 100 * factor)

stress = []
strain = []

plt.figure(figsize=(16, 9))

# Initial stress-strain calculation using the first element of the sequence
stress.append(sl.solve_Neuber(0, 0, sequence[0], Kt, E, Ku, nu, 'R-O'))
strain.append(sl.ramberg_osgood(stress[0], E, Ku, nu))
sl.ramberg_osgood_plot(stress[0], E, Ku, nu)

s_ro = np.linspace(0, np.max(sequence), 100 * factor)
e_ro = s_ro / E + (s_ro / Ku) ** (1 / nu)

s_mas = []
r_o_current = [0, 0]
loop = []

i = 1
smas = []
emas = []
while i < len(sequence):
    S = sequence[i]
    DS_temp = S - sequence[i - 1]
    DS = abs(DS_temp)

    Ds = sl.solve_Neuber(sequence[i - 1], strain[-1], S, Kt, E, Ku, nu, 'Masing')
    De = sl.masing(Ds, E, Ku, nu)



    direction = "down" if DS_temp < 0 else "up"

    if direction == "down":
        s = stress[i - 1] - Ds
        e = strain[i - 1] - De
        if s_mas:
            for k in range(len(s_mas)):
                if s < s_mas[k][0]:
                    if s <= min([mas[0] for mas in s_mas]):
                        for z in range(1, len(s_mas), 2):
                            # Ensure that s_mas[z] exists before accessing it
                            if z < len(s_mas):
                                sl.masing_plot(s_mas[z][2], s_mas[z][3], s_mas[z][0], E, Ku, nu, direction)
                                loop.append([s_mas[z][0], s_mas[z][1], s_mas[z][2], s_mas[z][3], (s_mas[z][0] + s_mas[z][2]) / 2, 1])

                        smas, emas = sl.masing_plot(s_mas[0][2], s_mas[0][3], s, E, Ku, nu, direction)
                        s_mas[0] = [smas[-1], emas[-1], s_mas[0][2], s_mas[0][3], 'direction']
                        del s_mas[1:]
                    else:
                        for z in range(k + 1, len(s_mas), 2):
                            # Ensure that s_mas[z] exists before accessing it
                            if z < len(s_mas):
                                sl.masing_plot(s_mas[z][2], s_mas[z][3], s_mas[z][0], E, Ku, nu, direction)
                                loop.append([s_mas[z][0], s_mas[z][1], s_mas[z][2], s_mas[z][3], (s_mas[z][0] + s_mas[z][2]) / 2,1])
                        smas, emas = sl.masing_plot(s_mas[k][2], s_mas[k][3], s, E, Ku, nu, direction)
                        s_mas[k] = [smas[-1], emas[-1], s_mas[k][2], s_mas[k][3], 'direction']
                        del s_mas[k + 1:]
                    break
                elif s == s_mas[k][0] and e == s_mas[k][1] and direction != s_mas[k][4] and stress[i - 1] == s_mas[k][
                    2] and strain[i - 1] == s_mas[k][3]:
                    smas, emas = sl.masing_plot(stress[i - 1], strain[i - 1], s, E, Ku, nu, direction)
                    loop.append([s, e, stress[i - 1], strain[i - 1], (s + stress[i - 1]) / 2, 1])
                    del s_mas[k]
                    break
                elif k == len(s_mas) - 1:
                    smas, emas = sl.masing_plot(stress[i - 1], strain[i - 1], s, E, Ku, nu, direction)
                    s_mas.append([smas[-1], emas[-1], smas[0], emas[0], direction])
        else:
            smas, emas = sl.masing_plot(stress[i - 1], strain[i - 1], s, E, Ku, nu, direction)
            s_mas.append([smas[-1], emas[-1], smas[0], emas[0], direction])
        s = smas[-1]
        e = emas[-1]

    elif direction == "up":
        s = stress[i - 1] + Ds
        e = strain[i - 1] + De

        if s >= r_o_current[0]:
            if s_mas:
                for k in range(len(s_mas)):
                    if s > s_mas[k][2] and e > s_mas[k][3] and direction != s_mas[k][4] and stress[i - 1] == s_mas[k][
                        0] and strain[i - 1] == s_mas[k][1]:
                        for z in range(2, len(s_mas)):
                            # Ensure that s_mas[z] exists before accessing it
                            if z < len(s_mas):
                                sl.masing_plot(s_mas[z][0], s_mas[z][1], s_mas[z][2], E, Ku, nu, direction)
                                loop.append(
                                    [s_mas[z][2], emas[-1], s_mas[z][0], s_mas[z][1], (s_mas[z][0] + s_mas[z][2]) / 2, 1])
                        loop.append([s_mas[0][0], s_mas[0][1], s_mas[0][2], s_mas[0][3], (s_mas[0][0] + s_mas[0][2]) / 2, 1])
                        smas, emas = sl.masing_plot(s_mas[0][0], s_mas[0][1], s_mas[0][2], E, Ku, nu, direction)
                        break
                    elif s == s_mas[k][2] and e == s_mas[k][3] and direction != s_mas[k][4] and stress[i - 1] == \
                            s_mas[k][0] and strain[i - 1] == s_mas[k][1]:
                        sl.masing_plot(stress[i - 1], strain[i - 1], s_mas[k][2], E, Ku, nu, direction)
                        smas, emas = sl.masing_plot(s_mas[k - 1][0], s_mas[k - 1][1], s_mas[0][2], E, Ku, nu, direction)
                        loop.append(
                            [smas[0], emas[0], smas[-1], emas[-1], (s_mas[k][0] + s_mas[k][2]) / 2, 1])

            s = sl.solve_Neuber(stress[i - 1], strain[i - 1], S, Kt, E, Ku, nu, 'R-O')
            sl.ramberg_osgood_plot(s, E, Ku, nu)
            e = sl.ramberg_osgood(s, E, Ku, nu)
            s_mas = []
            r_o_current = [s, e]
        else:
            for k in range(len(s_mas)):
                if s > s_mas[k][2]:
                    for z in range(k + 1, len(s_mas), 2):
                        # DEBUG
                        # Ensure that s_mas[z] exists before accessing it
                        if z < len(s_mas):
                            smas, emas = sl.masing_plot(s_mas[z][0], s_mas[z][1], s_mas[z][2], E, Ku, nu, direction)
                            loop.append(
                                [s_mas[z][2], emas[-1], s_mas[z][0], s_mas[z][1], (s_mas[z][0] + s_mas[z][2]) / 2,1])
                    smas, emas = sl.masing_plot(s_mas[k][0], s_mas[k][1], s, E, Ku, nu, direction)
                    del s_mas[k:]
                    s_mas.append([smas[0], emas[0], smas[-1], emas[-1], direction])
                    break
                elif s == s_mas[k][2] and e == s_mas[k][3] and direction != s_mas[k][4] and stress[i - 1] == s_mas[k][
                    0] and strain[i - 1] == s_mas[k][1]:
                    smas, emas = sl.masing_plot(stress[i - 1], strain[i - 1], s_mas[k][2], E, Ku, nu, direction)
                    loop.append([s_mas[k][0], s_mas[k][1], s_mas[k][2], s_mas[k][3], (s_mas[k][0] + s_mas[k][2]) / 2, 1])
                    del s_mas[k]
                    break
                elif k == len(s_mas) - 1:
                    smas, emas = sl.masing_plot(stress[i - 1], strain[i - 1], s, E, Ku, nu, direction)
                    s_mas.append([smas[0], emas[0], smas[-1], emas[-1], direction])
            s = smas[-1]
            e = emas[-1]
    stress.append(s)
    strain.append(e)
    i += 1





def solve_for_n(ea, sf, E, b, ef, c):
    def equation(N):
        return (sf / E) * (2 * N) ** b + ef * (2 * N) ** c - ea

    n_initial_guess = np.array(1e10)
    n_solution = opt.fsolve(equation, n_initial_guess, xtol=1e-10, maxfev=100000)

    return n_solution[0]


damage = 0
for i in range(len(loop)):
    ea = abs(loop[i][3] - loop[i][1])
    cycles_f = solve_for_n(ea, sf, E, b, ef, c)
    damage += 1 / cycles_f

print(f"Damage without mean stress effect: {(1 / damage):.3e}")


# me epidrash sm (PSWT)
def pswt_for_n(ea, smax, sf, E, b, ef, c):
    def equation_pswt(N):
        return (sf ** 2 / E) * (2 * N) ** (2 * b) + sf * ef * E * (2 * N) ** (b + c) - ea * E * smax

    n_initial_guess_pswt = np.array(1e8)
    n_solution_pswt = opt.fsolve(equation_pswt, n_initial_guess_pswt, xtol=1e-10, maxfev=100000)

    return n_solution_pswt[0]


damage_pswt = 0
for i in range(len(loop)):
    ea_pswt = abs(loop[i][3] - loop[i][1])
    smax = loop[i][4] + abs(loop[i][2] - loop[i][0])
    cycles_f_pswt = pswt_for_n(ea_pswt, smax, sf, E, b, ef, c)
    damage_pswt += 1 / cycles_f_pswt

print(f"Damage with mean stress effect - PSWT: {(1 / damage_pswt):.10e}")


# me epidrash sm (morrow)
def morr_for_n(de, ss, E, b, ef, c):
    def equation_morr(N):
        return (ss / E) * (2 * N) ** b + ef * (2 * N) ** c - de / 2

    n_initial_guess_morr = np.array(1e5)
    n_solution_morr = opt.fsolve(equation_morr, n_initial_guess_morr, xtol=1e-10, maxfev=100000)
    return n_solution_morr[0]


damage_morr = 0
for i in range(len(loop)):
    de = abs(loop[i][3] - loop[i][1])
    ss = sf-loop[i][4]
    cycles_f_morr = morr_for_n(de, ss, E, b, ef, c)
    damage_morr += 1 / cycles_f_morr
    # print(cycles_f_morr)

print(f"Damage with mean stress effect - MORROW: {(1 / damage_morr):.10e}")

plt.grid()
plt.show()
