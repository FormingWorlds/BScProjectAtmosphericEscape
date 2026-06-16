### This file will give the energy limited rate of escape for an atmosphere ###

import numpy as np
import scipy as sp
from rr_escape import calc_sound_speed, calc_sonic_point_radius

def el_escape_rate(epsilon_xuv, F_xuv, R_base, M_p, K_tide=1):
    '''
    This function returns the energy limited escape rate for an atmosphere following the formulation in Lopez 2017 (eq 2).


    Inputs:
    epsilon_xuv (float): species dependent evaporation efficiency [unitless]
    F_xuv (float): XUV flux [W/m^2]
    R_base (float): base radius of the atmosphere [m]
    M_p (float): planetary mass [kg]
    K_tide (float): tide constant [unitless]. 

    Output: 
    el_escape_rate (float): energy limited escape rate [kg/s]
    '''
    
    M_el_rate = epsilon_xuv * np.pi * F_xuv * R_base**3 / (sp.constants.G * M_p * K_tide)

    return M_el_rate

evaporation_efficiencies = {'H2': 0.1, #From Lopez 2017. Assumed to be 10% for Hydrogen and Helium, lighter elements, and 1% for heavier mixes
                            'H2He': 0.1,
                            'H2O':  0.01,
                            'CO2':  0.01,  
                            'N2': 0.01} 