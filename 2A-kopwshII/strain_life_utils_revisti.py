import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve


def solve_Neuber(DS,strain_from,stress_from,E,Ku,nu,behavior):

    if behavior == "R-O":
        DS = DS
        DE = DS/E
        e = strain_from + DE
        S = stress_from + DS
        es = S*e #*(Kt)**2

        eq = lambda s : es/(s+1e-10) - s/E - (s/Ku)**(1/nu)
        s_init = S

        s_final =fsolve(eq, s_init)
        
        # s_plot = np.linspace(stress_from, S, 100)
        # e_plot = es/s_plot
        # plt.plot(e_plot,s_plot)

        # DS_plot = np.linspace(0, DS, 100)
        # plt.plot(strain_from+DS_plot/E, stress_from+DS_plot)
        # print(stress_from,strain_from,S,e)
        # print(s_final)

        return s_final[0]

    elif behavior == "Masing":
        DeDs = (DS**2)/E

        eq = lambda Ds_r : DeDs/(Ds_r+1e-10) - Ds_r/E - 2*(Ds_r/(2*Ku))**(1/nu)
        Ds_init = DS

        Ds_final =fsolve(eq, Ds_init)

        # DS = DS*dir
        # DE = DS/E
        # e = strain_from + DE
        # S = stress_from + DS
        # DeDs = DS*DE #*(Kt)**2        
        # s_plot = np.linspace(stress_from+Ds_final[0]*dir, S, 100)
        # e_plot = strain_from + DeDs/(s_plot - stress_from)
        # plt.plot(e_plot,s_plot)

        # DS_plot = np.linspace(0, DS, 100)
        # plt.plot(strain_from+DS_plot/E, stress_from+DS_plot)

        return Ds_final[0]


def ramberg_osgood(s,E,Ku,nu):

    return s/E + (s/Ku)**(1/nu)

def ramberg_osgood_plot(stress,E,Ku,nu):

    s = np.linspace(0,stress,100)
    eps = ramberg_osgood(s,E,Ku,nu)
    # plt.yscale("log")
    # plt.xscale("log")
    plt.plot(eps,s)

def masing(Ds, E, Ku, nu):

    De_pl = 2*((Ds/2)/Ku)**(1/nu)
    De_el = Ds/E
    De = De_el + De_pl
    return De

def masing_plot(stress_from,strain_from,stress_to,E,Ku,nu,dir):
    s = np.linspace(stress_from,stress_to,100)

    Ds = abs(s - stress_from)
    De = masing(Ds,E,Ku,nu)

    e = strain_from + De * dir
    # plt.yscale("log")
    # plt.xscale("log")
    plt.plot(e,s)
