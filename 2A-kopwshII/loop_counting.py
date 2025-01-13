import strain_life_utils as sl
import numpy as np
from matplotlib import pyplot as plt
import utilities as utils


def main():


    E = 210000
    Ku = 1434
    nu=0.14
    Kt = 1.574
    s_normalization = 174
    factor = 100
    sequence = utils.reader("Strain-Gauge-Console")
    sequence = np.array(sequence[:20])*s_normalization/100 * Kt * factor #notch
    # print(sequence)

    # sequence = np.array([0.7, 0.3, 0.65, 0.2, 0.82, 0.5, 0.8,0.3, 1, 0.4])*s_normalization*factor
    err = 1/factor/100 # tune based on the factor
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
    s_ro = np.linspace(0,(np.max(sequence)),1000000)
    e_ro = s_ro/E + (s_ro/Ku)**(1/nu)
    # sel = [0, Kt*sequence[0]]
    # eel = [0,Kt*sequence[0]/E]
    # eneu = np.linspace(0.004,0.01,1000)
    # # print(sneu)
    # sneu = sel[-1] * eel[-1] /eneu
    # plt.plot(eneu,sneu)
    # plt.plot(eel, sel)
    # plt.plot(e_ro,s_ro)
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
            s_c = stress[i] + Ds
            e_c = strain[i] + De
        #
            index = 0  # Initialize the index
            for j in range(len(e_ro)):
                if abs(e_c - e_ro[j]) < err:  # Check if the absolute difference is small enough
                    index = j
                    break  # Stop searching once the first close value is found

            # Ensure an index was found
            if s_c >= s_ro[index]:  # Check if s_c is greater than or equal to s_ro at the found index
                s = s_ro[index]
                e = e_ro[index]
                sl.ramberg_osgood_plot(s,E,Ku,nu)
                sl.masing_plot(stress[i], strain[i], s-err2, E, Ku, nu, direction)

            else:
                s = s_c
                e = e_c
                sl.masing_plot(stress[i], strain[i], s, E, Ku, nu, direction)


        #
        stress.append(s)
        strain.append(e)


        # # branches.append(DS)
        # #
        # # if branches[i - 1] < DS:
        # #     # print('test',i)
        # #     DS = branches[i - 1]
        #
        # Ds = sl.solve_Neuber(DS,Kt,E,Ku,nu,'Masing')
        # De = sl.masing(Ds,E,Ku,nu)
        #

        # stress.append(s)
        # strain.append(9999)
        # sl.masing_plot(stress[i], strain[i], s, E, Ku, nu, direction)
        # s_check = s
        # e_check = e
        # for j in range(len(e_ro)):
        #     if abs(e_check-e_ro[j])<0.05:
        #         index = j
        #         break
        #
        # if s_check >= s_ro[index]:
        #     s = s_ro[index]
        #     stress.append(sl.solve_Neuber(S, Kt, E, Ku, nu))
        #     sl.ramberg_osgood_plot(stress[i],E,Ku, nu)
        #     strain.append(sl.ramberg_osgood(stress[i], E, Ku, nu))
        #     sl.masing_plot(stress[i],strain[i],s,E,Ku,nu,direction)

    # print(stress,strain)
    # print(np.array(branches))


    plt.grid()
    plt.show()
    # print(stress[:4])
if __name__ == "__main__":
    main()