import strain_life_utils as sl
import numpy as np
from matplotlib import pyplot as plt
import utilities as utils
import time
start_time = time.time()

# parameters
E = 210000
Ku = 1434
nu=0.14
Kt = 1.574
s_normalization = 174
factor = 1 # scale

# sequence read
sequence = utils.reader("Strain-Gauge-Console")
sequence = np.array(sequence[:30])*s_normalization/100 * factor #notch
# print(sequence)
# sequence = np.array([0.7, 0.3, 0.65, 0.2, 0.82, 0.5, 0.8,0.3, 1, 0.4])*s_normalization*factor

# def errors for loop counting
err_e = 1e-5 # tune based on maximum stress
err_s = 1e-2

stress = []
strain = []

plt.figure(figsize=(16, 9))

# initial R-O
stress.append(sl.solve_Neuber(0,0,sequence[0],Kt,E,Ku,nu,'R-O'))
strain.append(sl.ramberg_osgood(stress[0],E,Ku,nu))
# print(stress[0])
sl.ramberg_osgood_plot(stress[0],E,Ku,nu)
# sl.ramberg_osgood_plot(np.max(sequence), E,Ku,nu)

# initialize limit for Memory type 1
s_ro = np.linspace(0,(np.max(sequence)),100*factor)
e_ro = s_ro/E + (s_ro/Ku)**(1/nu)
s_c = 0
e_c = 0

s_mas = []
e_mas = []

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
        sl.masing_plot(stress[i], strain[i], s, E, Ku, nu, direction)
        smas, emas = sl.masing_curve(stress[i], strain[i], s, E, Ku, nu, direction)

    elif direction == "up":
        s = stress[i] + Ds
        e = strain[i] + De
        # print(stress[i], strain[i])
        smas, emas = sl.masing_curve(stress[i], strain[i], s, E, Ku, nu, direction)

        for smas_point, emas_point in zip(smas, emas):
            ro_index = np.argmin(np.abs(s_ro - smas_point))
            ro_strain = e_ro[ro_index]

            if emas_point <= ro_strain:
                Ds1 = smas_point - stress[i]
                DsN = Ds-Ds1
                # print(f"Intersection detected at stress: {smas_point}, strain: {emas_point}")
                plt.scatter(emas_point, smas_point, color='red', label='Intersection Point')
                s_c = smas_point
                e_c = emas_point

                s, e = sl.solve_Neuber_rec(smas[-1], smas_point, emas_point, Kt, E,Ku,nu)
                sl.ramberg_osgood_plot(s,E,Ku,nu)
                break

            else:
                s_c = s
                e_c = e


        sl.masing_plot(stress[i], strain[i], s_c, E, Ku, nu, direction)

    stress.append(s)
    strain.append(e)
    s_mas.append(smas)
    e_mas.append(emas)


# improved rainflow
loops = 0
loop = []
for j in range(1, len(s_mas)+1):

    min_s = np.min(s_mas[j - 1]) if e_mas[j - 1][0] - e_mas[j - 1][1] < 0 else np.max(s_mas[j - 1])
    max_s = np.max(s_mas[j - 1]) if e_mas[j - 1][0] - e_mas[j - 1][1] < 0 else np.min(s_mas[j - 1])
    min_e = np.min(e_mas[j - 1]) if e_mas[j - 1][0] - e_mas[j - 1][1] < 0 else np.max(e_mas[j - 1])
    max_e = np.max(e_mas[j - 1]) if e_mas[j - 1][0] - e_mas[j - 1][1] < 0 else np.min(e_mas[j - 1])

    # Combine into min and max pairs
    min_pair = np.array((min_s, min_e))
    max_pair = np.array((max_s, max_e))

    dir_excl = "d" if s_mas[j-1][0]-s_mas[j-1][1] > 0 else "u"


    for k in range(len(s_mas)):
        if k == j - 1:  # Skip the current (j-1)-th Masing curve
            continue

        # Combine stress and strain for the k-th Masing curve
        ch_temp = np.column_stack((s_mas[k], e_mas[k]))
        dir_temp = "d" if s_mas[k][0]-s_mas[k][1] > 0 else "u"


        # Check if (s, e)_min and (s, e)_max exist in the k-th Masing curve
        # Check if both components (s and e) of (s, e)_min are close to any pair in the k-th Masing curve
        ch1 = np.isclose(ch_temp[:, 0], min_pair[0], err_s) & np.isclose(ch_temp[:, 1], min_pair[1],err_e)

        # Check if both components (s and e) of (s, e)_max are close to any pair in the k-th Masing curve
        ch2 = np.isclose(ch_temp[:, 0], max_pair[0], err_s) & np.isclose(ch_temp[:, 1], max_pair[1],err_e)

        # print(ch1)
        # DEBUG: Print the comparison results for inspection
        # print(f"Checking curve {k} against curve {j - 1}")
        # print(dir_temp)
        # print(dir_excl)
        # print(f"Min Pair: {min_pair}, Max Pair: {max_pair}")
        # print(f"ch1: {np.any(ch1)}, ch2: {np.any(ch2)}")


        if np.any(ch1) == 1 and np.any(ch2) == 1:
            if dir_temp != dir_excl:
                if k > j-1:
                    # print(f"Curve {k} intersects both min and max of curve {j-1}")
                    # print(f"min pair: {min_pair}, max_pair: {max_pair}")
                    loops += 1
                    loop.append([j-1, min_pair[0], min_pair[1], max_pair[0], max_pair[1]])
                elif k < j-1:
                    if dir_temp == "d" and max_s < np.max(s_mas[k]):
                            # print(f"Curve {k} intersects both min and max of curve {j - 1}")
                            # print(f"min pair: {min_pair}, max_pair: {max_pair}")
                            loops += 1
                            loop.append([j - 1, min_pair[0], min_pair[1], max_pair[0], max_pair[1]])
                    else:
                        continue
                    if dir_temp == "u" and min_s > np.min(s_mas[k]):
                            # print(f"Curve {k} intersects both min and max of curve {j - 1}")
                            # print(f"min pair: {min_pair}, max_pair: {max_pair}")
                            loops += 1
                            loop.append([j - 1, min_pair[0], min_pair[1], max_pair[0], max_pair[1]])
                    else:
                        continue
            else:
                continue
# print(len(loop))

# REMOVE DOUBLE LOOPS
for z in range(len(loop) - 2, -1, -1):  # Iterate in reverse
    if (
        loop[z][0] == loop[z+1][0] - 1
        and loop[z][1] == loop[z+1][3]
        and loop[z][2] == loop[z+1][4]
        and loop[z][3] == loop[z+1][1]
        and loop[z][4] == loop[z+1][2]
    ):
        del loop[z+1]

# print(len(loop))












total_loops = len(loop)
print(total_loops)

plt.grid()
plt.show()

print(loop)












print("time elapsed: {:.2f}s".format(time.time() - start_time))