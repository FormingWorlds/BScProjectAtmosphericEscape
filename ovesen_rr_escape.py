'''
Malina Ovesen

Radiation-recombination-limited escape mechanism.
'''

import numpy as np
import scipy as sp
from ovesen_conversions import *
from ovesen_atmospheres.simple_H2 import T_wind, mu_wind, M_p, F_xuv, nu_0, mu_plus_wind, R_base

def find_R_base():
    #placeholder for finding R_base based on input parameters
    R_base = None
    return R_base


#In the input box I put the so-far necessary inputs from atmosphere profile for this func to work
def rr_escape_rate(T_wind, mu_wind, M_p, F_xuv, nu_0, mu_plus_wind, R_base):
    '''
    Calculates the radiation-recombination-limited escape rate.

    Takes input parameters: x [unit], y [unit], z [unit], ...

    All calculations done in SI units.
    '''
    #What we want to calculate and return
    #M_rr_rate = -4 * pi* rho_s * c_s * R_s**2   REFERENCE: Lopez 2017 eq. 4
    
    #REFERENCE: Murray et. al. 2009 after eq. 7
    alpha_rec_B = 2.7 * 10**(-13) * (T_wind / 10**4)**0.9 * cm_to_m**3 #[m^3 s^-1] recombination coefficient for case B recombination, where T_wind: temperature of the escaping atmosphere

    ### calculates c_s based on input ###
    k_b = sp.constants.k #[J K^-1]
    m_p = sp.constants.m_p #[kg]
    c_s = np.sqrt((k_b * T_wind) / (mu_wind * m_p)) #[m s^-1] sound speed, where k_B: Boltzmann constant, T: temperature, mu: mean molecular weight, m_p: proton mass
    

    ### calculates rho_s based on input ###
    G = sp.constants.G #[m^3 kg^-1 s^-2] 
    R_s = G * M_p / (2 * c_s**2) #[m] radius to the sonic point, where G: gravitational constant, M_p: planetary mass, c_s: sound speed
    #!!!need to remember to check if R_s > R_base, where R_base is the radius of the base of the escaping atmosphere, otherwise the escape rate is not radiation-recombination-limited (?). I can set R_s = R_base if smaller
    

    ### finds R_base based on input ###
    #placeholder


    ### calculates rho_s based on input ### 
    h = sp.constants.h #[J s] 
    # reference formula: ovesen math calc derivation. Needed to place n_0_base straight into n_plus_base to cancel sigma_nu0
    n_plus_base = np.sqrt( F_xuv * G * M_p / (h * nu_0 * alpha_rec_B * c_s**2 * R_base**2) ) #[m^-3] number density of the escaping atmosphere at the base, where F_xuv: XUV flux, G: gravitational constant, M_p: planetary mass, h: Planck's constant, nu_0: frequency of the ionising radiation, alpha_rec_B: recombination coefficient for case B recombination, c_s: sound speed, R_base: radius of the base of the escaping atmosphere

    rho_base = n_plus_base * mu_plus_wind * m_p #[kg m^-3] density at the base of the escaping atmosphere, where n_plus_base: number density of the escaping atmosphere at the base, mu_wind: mean molecular weight, m_p: proton mass
    
    #REFERENCE Lopez 2017 eq. 5
    rho_s = rho_base * np.exp((-G * M_p) / (R_base * c_s**2) * (R_base/R_s - 1 )) #[kg m^-3] rho_s is the density at the sonic point, rho_base: density at the base of the escaping atmosphere, R_base: radius of the base of the escaping atmosphere

    
    M_rr_rate = -4 * np.pi* rho_s * c_s * R_s**2 #[kg s^-1] radiation-recombination-limited escape rate, where rho_s: density at the sonic point, c_s: sound speed, R_s: radius to the sonic point

    return M_rr_rate 

print(rr_escape_rate(T_wind, mu_wind, M_p, F_xuv, nu_0, mu_plus_wind, R_base))