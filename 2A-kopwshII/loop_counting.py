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

    sequence = np.array([ 25, 5, 14, -14, 30])*20

    stress = []
    strain = []
    branches = []
    plt.figure(figsize=(20, 6))
    stress.append(sl.solve_Neuber(sequence[0],Kt,E,Ku,nu,'R-O'))
    strain.append(sl.ramberg_osgood(stress[0],E,Ku,nu))

    sl.ramberg_osgood_plot(stress[0],E,Ku,nu)

    for i,S in enumerate(sequence[1:]):

        DS_temp = S - sequence[i]
        DS = abs(DS_temp)

        if DS_temp < 0:
            direction = "down"
        elif DS_temp > 0:
            direction = "up"

        branches.append(DS)

        if branches[i-1] < DS:
            # print('test',i)
            DS = branches[i-1]


        Ds = sl.solve_Neuber(DS,Kt,E,Ku,nu,'Masing')
        De = sl.masing(Ds,E,Ku,nu)

        if direction == "down":
            s = stress[i] - Ds
            e = strain[i] - De
        elif direction == "up":
            s = stress[i] + Ds
            e = strain[i] + De


        # if s == stress[i]:


        stress.append(s)
        strain.append(e)
        sl.masing_plot(stress[i],strain[i],s,E,Ku,nu,direction)

    print(stress,strain)
    print(branches)



    plt.show()
if __name__ == "__main__":
    main()