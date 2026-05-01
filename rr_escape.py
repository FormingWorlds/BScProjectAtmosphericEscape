'''
Malina Ovesen

Radiation-recombination-limited escape mechanism.
'''

import numpy as np
import scipy as sp
from conversions import *

def find_R_base(P_base, radii, pressures):
    '''
    Finds the radius of the base of the escaping atmosphere.

    Takes input parameters: x [unit], y [unit], z [unit], ...

    All calculations done in SI units.
    '''
    #finds the absolute difference between the pressure profile and the target pressure at the base of the escaping atmosphere
    difference_array = np.absolute(pressures - P_base) #[Pa] array of the absolute difference between the pressure profile and the pressure at the base of the escaping atmosphere, where pressures: pressure profile of the atmosphere based on the barometric formula, P_base: pressure at the base of the escaping atmosphere

    #finds the index of minimum element from the array
    index = difference_array.argmin()
    R_base = radii[index] #[m] radius of the base of the escaping atmosphere, where radii: array of radii from the planetary radius to 10 times the planetary radius, index: index of minimum element from the array of the absolute difference between the pressure profile and the pressure at the base of the escaping atmosphere
  
    return R_base


def calc_sound_speed(T_wind, mu_wind):
    '''
    Calculates the sound speed of the escaping atmosphere.

    Takes input parameters: x [unit], y [unit], z [unit], ...

    All calculations done in SI units.
    '''
    k_b = sp.constants.k #[J K^-1]
    m_p = sp.constants.m_p #[kg]
    c_s = np.sqrt((k_b * T_wind) / (mu_wind * m_p)) #[m s^-1] sound speed, where k_B: Boltzmann constant, T: temperature, mu: mean molecular weight, m_p: proton mass

    return c_s


def calc_sonic_point_radius(M_p, c_s, R_base):
    '''
    Calculates the radius to the sonic point.

    Takes input parameters: x [unit], y [unit], z [unit], ...

    All calculations done in SI units.
    '''
    ### calculates R_s based on input ###
    G = sp.constants.G #[m^3 kg^-1 s^-2] 
    R_s = G * M_p / (2 * c_s**2) #[m] radius to the sonic point, where G: gravitational constant, M_p: planetary mass, c_s: sound speed
    print("Calculated R_s:", R_s, "m")
    # checks if R_s is smaller than R_base, if so, sets R_s = R_base and prints a message, otherwise keeps R_s = G * M_p / (2 * c_s**2) and prints a message
    if R_s < R_base:
        print("R_s is smaller than R_base, escape is not radiation-recombination-limited. Setting R_s = R_base.")
        R_s = R_base
    else: 
        print("R_s is larger than R_base, escape is radiation-recombination-limited. Keeping R_s = G * M_p / (2 * c_s**2).")
        
    return R_s


def calc_density_at_sonic_point(M_p, F_xuv, nu_0, T_wind, R_s, c_s, R_base, mu_plus_wind):
    '''
    Calculates the density at the sonic point.

    Takes input parameters: x [unit], y [unit], z [unit], ...

    All calculations done in SI units.
    '''

    h = sp.constants.h #[J s] 
    #REFERENCE: Murray et. al. 2009 after eq. 7
    alpha_rec_B = 2.7 * 10**(-13) * (T_wind / 10**4)**0.9 * cm_to_m**3 #[m^3 s^-1] recombination coefficient for case B recombination, where T_wind: temperature of the escaping atmosphere

    G = sp.constants.G #[m^3 kg^-1 s^-2] 
    # reference formula: ovesen math calc derivation. Needed to place n_0_base straight into n_plus_base to cancel sigma_nu0
    n_plus_base = np.sqrt( F_xuv * G * M_p / (h * nu_0 * alpha_rec_B * c_s**2 * R_base**2) ) #[m^-3] number density of the escaping atmosphere at the base, where F_xuv: XUV flux, G: gravitational constant, M_p: planetary mass, h: Planck's constant, nu_0: frequency of the ionising radiation, alpha_rec_B: recombination coefficient for case B recombination, c_s: sound speed, R_base: radius of the base of the escaping atmosphere

    m_p = sp.constants.m_p #[kg]
    rho_base = n_plus_base * mu_plus_wind * m_p #[kg m^-3] density at the base of the escaping atmosphere, where n_plus_base: number density of the escaping atmosphere at the base, mu_wind: mean molecular weight, m_p: proton mass
    
    #REFERENCE Lopez 2017 eq. 5
    rho_s = rho_base * np.exp((G * M_p) / (R_base * c_s**2) * (R_base/R_s - 1 ) - 1/2) #[kg m^-3] rho_s is the density at the sonic point, rho_base: density at the base of the escaping atmosphere, R_base: radius of the base of the escaping atmosphere

    return rho_s

#In the input box I put the so-far necessary inputs from atmosphere profile for this func to work
def rr_escape_rate(rho_s, c_s, R_s):
    '''
    Calculates the radiation-recombination-limited escape rate.

    Takes input parameters: x [unit], y [unit], z [unit], ...

    All calculations done in SI units.
    '''
 
    M_rr_rate = -4 * np.pi* rho_s * c_s * R_s**2 #[kg s^-1] radiation-recombination-limited escape rate REFERENCE: Lopez 2017 eq. 4, where rho_s: density at the sonic point, c_s: sound speed, R_s: radius to the sonic point

    return M_rr_rate 

