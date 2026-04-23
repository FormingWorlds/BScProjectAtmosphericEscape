'''
Malina Ovesen

Radiation-recombination-limited escape mechanism.
'''

import numpy as np
import scipy as sp

#In the input box I put the so-far necessary inputs from atmosphere profile for this func to work
def rr_escape_rate(T_wind, mu_wind, M_p, F_xuv, sigma_nu0, nu_0, alpha_rec_B, mu_plus_wind, R_base):
    '''
    Calculates the radiation-recombination-limited escape rate.

    Takes input parameters: x [unit], y [unit], z [unit], ...

    All calculations done in SI units.
    '''
    #What we want to calculate and return
    #M_rr_rate = -4 * pi* rho_s * c_s * R_s**2   REFERENCE: Lopez 2017 eq. 4


    ### calculates c_s based on input ###
    k_b = sp.constants.k 
    m_p = sp.constants.m_p
    c_s = np.sqrt((k_b * T_wind) / (mu_wind * m_p)) #sound speed, where k_B: Boltzmann constant, T: temperature, mu: mean molecular weight, m_p: proton mass
    

    ### calculates rho_s based on input ###
    G = sp.constants.G
    R_s = G * M_p / (2 * c_s**2) #radius to the sonic point, where G: gravitational constant, M_p: planetary mass, c_s: sound speed
    #!!!need to remember to check if R_s > R_base, where R_base is the radius of the base of the escaping atmosphere, otherwise the escape rate is not radiation-recombination-limited (?). I can set R_s = R_base if smaller
    

    ### calculates rho_s based on input ### 
    # reference formula: rho_s = rho_base * exp((-G * M_p) / (R_base * c_s**2) * (R_base/R_s - 1 )) 
    g_base = G * M_p / R_base**2 #gravitational acceleration at the base of the escaping atmosphere, where G: gravitational constant, M_p: planetary mass, R_base: radius of the base of the escaping atmosphere
    H_base = c_s**2 / g_base #scale height at the base of the escaping atmosphere, where c_s: sound speed, g_base: gravitational acceleration at the base of the escaping atmosphere

    n_0_base = 1 / (sigma_nu0 * H_base) #number density at the base of the escaping atmosphere, where sigma_nu0: photoionization cross section at the ionization threshold (?), H_base: scale height at the base of the escaping atmosphere

    h = sp.constants.h
    n_plus_base = np.sqrt( F_xuv * sigma_nu0 * n_0_base / (h * nu_0 * alpha_rec_B) )

    rho_base = n_plus_base * mu_plus_wind * m_p #density at the base of the escaping atmosphere, where n_plus_base: number density of the escaping atmosphere at the base, mu_wind: mean molecular weight, m_p: proton mass
    rho_s = rho_base * np.exp((-G * M_p) / (R_base * c_s**2) * (R_base/R_s - 1 )) #rho_s is the density at the sonic point, rho_base: density at the base of the escaping atmosphere, R_base: radius of the base of the escaping atmosphere

    return None # Placeholder for the actual implementation

