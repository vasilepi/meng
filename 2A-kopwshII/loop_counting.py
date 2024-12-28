from matplotlib.pyplot import figure

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


    # sequence = utils.reader("Strain-Gauge-Console")
    # sequence = np.array(sequence)*s_normalization/100 * Kt #notch
    # print(sequence)

    sequence = np.array([25, 5, 14, -14, 30, 10])*20

    stress = []
    strain = []
    branches = []
    plt.figure(figsize=(10, 6))
    stress.append(sl.solve_Neuber(0,0,sequence[0],Kt,E,Ku,nu,'R-O'))
    strain.append(sl.ramberg_osgood(stress[0],E,Ku,nu))
    # print(stress[0])
    sl.ramberg_osgood_plot(stress[0],E,Ku,nu)

    s_ro = np.linspace(0,(np.max(sequence)),1000)
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
        print(i)
        Ds = sl.solve_Neuber(sequence[i],strain[i],S,Kt,E,Ku,nu,'Masing')
        De = sl.masing(Ds,E,Ku,nu)
        #
        # if direction == "down":
        #     s = stress[i] - Ds
        #     e = strain[i] - De
        # elif direction == "up":
        #     s = stress[i] + Ds
        #     e = strain[i] + De
        #
        #     index = 0
        #     for j in range(len(e_ro)):
        #         if abs(e-e_ro[j]) < 5e-8:
        #             index = j
        #             break
        #
        #     if s >= s_ro[index]:
        #         s = sl.solve_Neuber(sequence[i],strain[i],S,Kt,E,Ku,nu,'R-O')
        #         e = sl.ramberg_osgood_plot(s,E,Ku,nu)
        #
        # stress.append(s)
        # strain.append(e)


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



    plt.show()
if __name__ == "__main__":
    main()