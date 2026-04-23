'''
Malina Ovesen

Radition-recombination-limited escape mechanism.
'''

import numpy as np
import scipy as sp

def rr_escape_rate():
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
    
    return None # Placeholder for the actual implementation

