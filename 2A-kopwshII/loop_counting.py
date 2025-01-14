import strain_life_utils as sl
import numpy as np
from matplotlib import pyplot as plt
import utilities as utils





E = 210000
Ku = 1434
nu=0.14
Kt = 1.574
s_normalization = 174
factor = 100
sequence = utils.reader("Strain-Gauge-Console")
sequence = np.array(sequence[:10])*s_normalization/100 * Kt * factor #notch
# print(sequence)

# sequence = np.array([0.7, 0.3, 0.65, 0.2, 0.82, 0.5, 0.8,0.3, 1, 0.4])*s_normalization*factor
err = 1/factor/1000 # tune based on the factor
err2 = s_normalization/(factor*1e-1) * Kt
stress = []
strain = []
branches = []
plt.figure(figsize=(10, 6))
stress.append(sl.solve_Neuber(0,0,sequence[0],Kt,E,Ku,nu,'R-O'))
strain.append(sl.ramberg_osgood(stress[0],E,Ku,nu))
# print(stress[0])
sl.ramberg_osgood_plot(stress[0],E,Ku,nu)
# sl.ramberg_osgood_plot(np.max(sequence), E,Ku,nu)
s_ro = np.linspace(0,(np.max(sequence)),100*factor)
e_ro = s_ro/E + (s_ro/Ku)**(1/nu)
s_c = 0
e_c = 0
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

plt.grid()
plt.show()
print(stress[-1], strain[-1], stress[-2], strain[-2])
