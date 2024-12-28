import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

    

def solve_Neuber(stress_from,strain_from,S,Kt,E,Ku,nu,behavior):

    if behavior == "R-O":
        DS = abs(S - stress_from)
        DE = DS/E
        e = strain_from + DE
        es = S*e*(Kt)**2

        eq = lambda s : es/(s+1e-10) - s/E - (s/Ku)**(1/nu)
        s_init = S

        s_final =fsolve(eq, s_init)
        return s_final[0]

    elif behavior == "Masing":
        DS = abs(S - stress_from)
        DeDs = ((Kt*DS)**2)/E

        eq = lambda Ds_r : DeDs/(Ds_r+1e-10) - Ds_r/E - 2*(Ds_r/(2*Ku))**(1/nu)
        Ds_init = DS

        Ds_final =fsolve(eq, Ds_init)

        return Ds_final[0]
        

def ramberg_osgood(s,E,Ku,nu):
    
    return s/E + (s/Ku)**(1/nu)

def ramberg_osgood_plot(stress,E,Ku,nu):
    
    s = np.linspace(0,stress,100)   
    eps = ramberg_osgood(s,E,Ku,nu)  
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
    
    if dir == "up":
        e = strain_from + De
    
    elif dir == "down":
        e = strain_from - De


    plt.plot(e,s)
