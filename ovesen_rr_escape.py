'''
Malina Ovesen

Radition-recombination-limited escape mechanism.
'''

import numpy as np
import scipy as sp

#In the input box I put the so-far necessary inputs from atmosphere profile for this func to work
def rr_escape_rate(T_wind, mu_wind, M_p):
    '''
    Calculates the radiation-recombination-limited escape rate.

    Takes input parameters: x [unit], y [unit], z [unit], ...

    All calculations done in SI units.
    '''

    #M_rr_rate = -4 * pi* rho_s * c_s * R_s**2   REFERENCE: Lopez 2017 eq. 4

    ### calculates c_s based on input ###
    k_b = sp.constants.k 
    m_p = sp.constants.m_p
    c_s = np.sqrt((k_b * T_wind) / (mu_wind * m_p)) #sound speed, where k_B: Boltzmann constant, T: temperature, mu: mean molecular weight, m_p: proton mass
    
    ### calculates rho_s based on input ###
    G = sp.constants.G
    R_s = G * M_p / (2 * c_s**2) #radius to the sonic point, where G: gravitational constant, M_p: planetary mass, c_s: sound speed
    
    return None # Placeholder for the actual implementation

