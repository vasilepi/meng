import strain_life_utils as sl
import numpy as np
from matplotlib import pyplot as plt
import utilities as utils
import time
import scipy.optimize as opt
from scipy.optimize import root_scalar
from scipy.optimize import fsolve
start_time = time.time()

E = 210000
Ku = 922
nu = 0.15
b = -0.0971
c = -0.619
sf = 774
ef = 0.502
Kt = 2.17
s_normalization = 174
factor = 100
sequence = utils.reader("Strain-Gauge-Console")
sequence = (np.array(sequence[:9]))*(s_normalization/100  * factor) #notch
# print(sequence)

# sequence = np.array([0.7, 0.3, 0.65, 0.2, 0.82, 0.5, 0.8,0.3, 1, 0.4])*s_normalization*factor
err_e = 1e-4 * factor # tune based on the factor
err_s = 1e-4 * factor
err3 = 2e-5 * factor
stress = []
strain = []

# initiate stress-strain graph
plt.figure(figsize=(16, 9))

stress.append(sl.solve_Neuber(0,0,sequence[0],Kt,E,Ku,nu,'R-O'))
strain.append(sl.ramberg_osgood(stress[0],E,Ku,nu))
# print(stress[0])
sl.ramberg_osgood_plot(stress[0],E,Ku,nu)
# sl.ramberg_osgood_plot(np.max(sequence), E,Ku,nu)
s_ro = np.linspace(0,(np.max(sequence)),100*factor)
e_ro = s_ro/E + (s_ro/Ku)**(1/nu)


s_mas = []
e_mas = []

loop = []
for i,S in enumerate(sequence[1:]):
    DS_temp = S - sequence[i]
    DS = abs(DS_temp)

    if DS_temp < 0:
        direction = "down"
    elif DS_temp > 0:
        direction = "up"
    # print(i)
    Ds = sl.solve_Neuber(sequence[i],strain[i],S,Kt,E,Ku,nu,'Masing')
    De = sl.masing(Ds,E,Ku,nu)
    if direction == "down":
        s = stress[i] - Ds
        e = strain[i] - De
        smas, emas = sl.masing_curve(stress[i], strain[i], s, E, Ku, nu, direction)
        # sl.masing_plot(stress[i], strain[i], s, E, Ku, nu, direction)
        # Check for intersections with previous Masing curves of the same direction
        # Intersection Handling
        d_intersected = False  # Track if an intersection occurs
        for k in range(len(s_mas)):
            if k >= i:  # Skip future curves
                continue

            # Match direction
            prev_direction = "down" if s_mas[k][0] - s_mas[k][1] > 0 else "up"
            if prev_direction != direction:
                continue

            # Check intersection
            for smas_point, emas_point in zip(smas, emas):
                prev_index = -1
                prev_strain = e_mas[k][prev_index]

                if np.isclose(emas_point, prev_strain, err_e) and np.isclose(smas_point, s_mas[k][prev_index],
                                                                                  err_s):
                    d_intersected = True
                    plt.scatter(emas_point, smas_point, color='blue', label='Masing Intersection Point')

                    # Update Neuber and recalculate curve
                    Ds_new = sl.solve_Neuber(sequence[k], e_mas[k][0], S, Kt, E, Ku, nu, 'Masing')
                    De_new = sl.masing(Ds_new, E, Ku, nu)
                    s = stress[k] - Ds_new
                    e = strain[k] - De_new

                    # Plot updated curve
                    sl.masing_plot(stress[i], strain[i], smas_point, E, Ku, nu, direction)
                    sl.masing_plot(stress[k], strain[k], s, E, Ku, nu, direction)
                    smas_oldnew, emas_oldnew = sl.masing_curve(stress[k], strain[k], s, E, Ku, nu, direction)
                    s_mas[k] = smas_oldnew
                    e_mas[k] = emas_oldnew
                    break  # Exit inner loop once intersection is handled
            if d_intersected:
                loop.append([smas_point, emas_point, stress[i], strain[i]])
                break  # Exit outer loop if intersection is handled

        # If no intersection occurred, plot original curve
        if not d_intersected:
            sl.masing_plot(stress[i], strain[i], s, E, Ku, nu, direction)



    elif direction == "up":
        s = stress[i] + Ds
        e = strain[i] + De
        smas, emas = sl.masing_curve(stress[i], strain[i], s, E, Ku, nu, direction)
        u_intersected = False

        # Handle Up Intersections (with Ramberg-Osgood and other Masing curves)


        # Check intersections with other "up" Masing curves
        for k in range(len(s_mas)):
            if k >= i:
                continue
            prev_direction_u = "down" if s_mas[k][0] - s_mas[k][1] > 0 else "up"
            if prev_direction_u != direction:
                continue

            for smas_point, emas_point in zip(smas, emas):
                prev_index = -1
                prev_strain = e_mas[k][prev_index]

                if np.isclose(emas_point, prev_strain, err_e) and np.isclose(smas_point, s_mas[k][prev_index],
                                                                             err_s):
                    u_intersected = True
                    plt.scatter(emas_point, smas_point, color='green', label='Up Masing Intersection Point')

                    Ds_new = sl.solve_Neuber(sequence[k], e_mas[k][0], S, Kt, E, Ku, nu, 'Masing')
                    De_new = sl.masing(Ds_new, E, Ku, nu)
                    s = stress[k] + Ds_new
                    e = strain[k] + De_new

                    sl.masing_plot(stress[i], strain[i], smas_point, E, Ku, nu, direction)


                    if s < np.max(stress):
                        sl.masing_plot(stress[k], strain[k], s, E, Ku, nu, direction)
                        smas_oldnew, emas_oldnew = sl.masing_curve(stress[k], strain[k], s, E, Ku, nu, direction)
                        s_mas[k] = smas_oldnew
                        e_mas[k] = emas_oldnew
                    elif s >= np.max(stress):
                        sl.masing_plot(stress[k], strain[k], np.max(stress), E, Ku, nu, direction)
                        smas_oldnew, emas_oldnew = sl.masing_curve(stress[k], strain[k], np.max(stress), E, Ku, nu, direction)
                        s_mas[k] = smas_oldnew
                        e_mas[k] = emas_oldnew
                        s = sl.solve_Neuber(sequence[k], e_mas[k][0], S, Kt, E, Ku, nu, 'R-O')
                        e = sl.ramberg_osgood(s, E, Ku, nu)
                        sl.ramberg_osgood_plot(s, E, Ku, nu)
                        loop.append([stress[k], strain[k], np.max(stress), np.max(strain)])
                    break
            if u_intersected:
                loop.append([stress[i], strain[i], smas_point, emas_point])
                break

        if not u_intersected:
            for smas_point, emas_point in zip(smas, emas):
                ro_index = np.argmin(np.abs(s_ro - smas_point))
                ro_strain = e_ro[ro_index]

                if emas_point - ro_strain <= err3:

                    plt.scatter(emas_point, smas_point, color='red', label='Ramberg-Osgood Intersection Point')
                    s_c = smas_point
                    e_c = emas_point

                    s, e = sl.solve_Neuber_rec(smas[-1], smas_point, emas_point, Kt, E, Ku, nu)
                    sl.ramberg_osgood_plot(s, E, Ku, nu)
                    loop.append([stress[i], strain[i], s_c, e_c])
                    break
                else:
                    s_c = s
                    e_c = e
            sl.masing_plot(stress[i], strain[i], s_c, E, Ku, nu, direction)

    stress.append(s)
    strain.append(e)
    s_mas.append(smas)
    e_mas.append(emas)


# cumulative count
count = {}
for arr in loop:
    arr_tuple = tuple(arr)
    if arr_tuple in count:
        count[arr_tuple] += 1
    else:
        count[arr_tuple] = 1

# Convert back to the desired output format
loops_cum = [[count, arr] for arr, count in count.items()]
# print(loops_cum)





# stress-strain
plt.grid()
plt.show()



# print(loop)
print(len(loop))

###### CYCLES TO FAILURE ######
# xwris epidrash sm
def solve_for_N(ea, sf, E, b, ef, c):
    def equation(N):
        return (sf / E) * (2 * N) ** b + ef * (2 * N) ** c - ea

    N_initial_guess = 1e10
    N_solution = opt.fsolve(equation, N_initial_guess, xtol=1e-15, maxfev=100000)

    return N_solution[0]


damage = 0
for i in range(len(loops_cum)):
    ea = abs(loops_cum[i][1][3] - loops_cum[i][1][1])
    cycles_f = solve_for_N(ea, sf, E, b, ef, c)
    damage += 1 / cycles_f
    
print(f"Damage without mean stress effect: {(1 / damage):.3e}")


# me epidrash sm (PSWT)
def PSWT_for_N(ea, smax, sf, E, b, ef, c):
    def equation_PSWT(N):
        return (sf**2 / E) * (2 * N) ** (2*b) + sf*ef*E * (2 * N) ** (b+c) - ea*E*smax

    N_initial_guess_pswt = 1e8
    N_solution_pswt = opt.fsolve(equation_PSWT, N_initial_guess_pswt, xtol=1e-15, maxfev=100000)

    return N_solution_pswt[0]


damage_pswt = 0
for i in range(len(loops_cum)):
    ea_pswt = abs(loops_cum[i][1][3] - loops_cum[i][1][1])
    smax = (loops_cum[i][1][2]+loops_cum[i][1][0])/2 + abs(loops_cum[i][1][2] - loops_cum[i][1][0])
    cycles_f_pswt = PSWT_for_N(ea_pswt, smax, sf, E, b, ef, c)
    damage_pswt += 1 / cycles_f_pswt

print(f"Damage with mean stress effect - PSWT: {(1 / damage_pswt):.3e}")


# me epidrash sm (morrow)
def morr_for_N(de, sm, sf, E, b, ef, c):
    def equation_morr(N):
        return ((sf-sm) / E) * (2 * N) ** b + ef * (2 * N) ** c - de/2

    N_initial_guess_morr = 1e8
    N_solution_morr = opt.fsolve(equation_morr, N_initial_guess_morr, maxfev=100000)

    return N_solution_morr[0]


damage_morr = 0
for i in range(len(loops_cum)):
    de = abs(loops_cum[i][1][3] - loops_cum[i][1][1])
    sm = (loops_cum[i][1][2]+loops_cum[i][1][0])/2
    cycles_f_morr = morr_for_N(de, sm, sf, E, b, ef, c)
    damage_morr += 1 / cycles_f_morr

print(f"Damage with mean stress effect - MORROW: {(1 / damage_morr):.3e}")





# strain-life
Nplot = np.arange(1e3,1e6, 1)
eaplot = (sf/E * (2*Nplot)**b + ef*(2*Nplot)**c)*100

# strain-life graph
# plt.figure(figsize = (10,6))
# plt.plot(Nplot, eaplot)
# plt.grid()
# plt.xscale('log')
# plt.ylim((0, 1))
# plt.show()

print("time elapsed: {:.2f}s".format(time.time() - start_time))