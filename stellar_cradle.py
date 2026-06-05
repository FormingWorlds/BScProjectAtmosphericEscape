### This file will give birth to different stellar types which will be used to investigate the effect of different stellar types on the escape regime and rate of an atmosphere. ###
import numpy as np
import scipy as sp
from matplotlib import pyplot as plt

solar_mass = 1.989 * 10**30 #kg
solar_radius = 6.957 * 10**8 #m
solar_luminosity = 3.828 * 10**26 #W

def mass_to_luminosity(mass):
    '''
    This function takes in the mass of a star and returns the luminosity of the star in using the mass-luminosity relation from Eker et. al. 2018.
    Their luminosity relation consits of 6 piecewise equations which tie different ranges of stellar masses to a luminosity. As the paper states, the relation is only valid for stars between 0.179 and 31 solar masses, '
    therefore this function only covers that range as well.


    Input: stellar mass [kg]
    Output: stellar luminosity [W]
    '''

    #Relation uses solar masses so must convert before using the relation
    
    mass_sm = mass / solar_mass

    if 0.179 < mass_sm <= 0.45:
        log_L = 2.028 * np.log10(mass_sm) - 0.976
    elif 0.45 < mass_sm <= 0.72:
        log_L = 4.572 * np.log10(mass_sm) - 0.102
    elif 0.72 < mass_sm <= 1.05:
        log_L = 5.743 * np.log10(mass_sm) - 0.007
    elif 1.05 < mass_sm <= 2.4:
        log_L = 4.329 * np.log10(mass_sm) + 0.010
    elif 2.4 < mass_sm <= 7:
        log_L = 3.967 * np.log10(mass_sm) + 0.093
    elif 7 < mass_sm <= 31:
        log_L = 2.865 * np.log10(mass_sm) + 1.105
    else:
        raise ValueError('Mass must be between 0.179 and 31 solar masses. Other masses are not covered by the relation from Eker et. al. 2018.')

    L = 10**log_L #[W]

    return L

def mass_to_radius(mass):
    '''
    This function takes the mass of a star and converts it to a radius using the mass-radius relation from Giménez and Zamorano (1985).
    
    Input: stellar mass [kg]
    Output: stellar radius [m]
    '''
    mass_sm = mass / solar_mass #[solar masses]
    
    if mass_sm < 1.8:
        radius_sm = 1.13 * mass_sm**0.98 #[solar radii]
    elif mass_sm >= 1.8:
        radius_sm = 1.42 * mass_sm**0.56 #[solar radii]
    else:
        raise ValueError('Mass must be a positive value.')
    
    radius = radius_sm * solar_radius #[m]

    return radius

def calculate_T_eff(L, R):
    '''
    This function takes the luminosity and radius of a star and calculates the effective temperature of the star using the Stefan-Boltzmann law.

    Input: luminosity [W], radius [m]
    Output: effective temperature [K]
    '''
    sigma = sp.constants.sigma #[W m^-2 K^-4] Stefan-Boltzmann constant

    T_eff = (L / (4 * np.pi * R**2 * sigma))**(1/4) #[K]

    return T_eff

########################################################################################################################
### Values from Table 6 for calibration of B-V colour for population 1 main sequence stars from Böhm-Vitense (1970) ###

# Radiative atmosphere
T_eff_radiative_array = np.array([9380, 9210, 9080, 8930, 8790, 8630, 8470, 8310, 8150, 7990, 7830, 7680, 
                                  7530, 7380, 7250, 7120, 7010, 6890, 6790, 6680, 6580, 6470, 6390, 6270])[::-1] #[K] array to store the effective temperatures of the stars
B_V_radiative_array = np.array([0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 
                                0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46])[::-1] #[mag] array to store the B-V colours of the stars


# Convective atmosphere
T_eff_convective_array = np.array([7870, 7730, 7600, 7480, 7350, 7230, 7120, 7010, 6910, 6800, 6700, 6590, 6480, 6380, 6280, 
                                   6180, 6080, 5980, 5890, 5800, 5730, 5650, 5570, 5500, 5390, 5230, 5070, 4930, 4790, 4660])[::-1] #[K] array to store the effective temperatures of the stars
B_V_convective_array = np.array([0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.50, 0.52, 0.54,
                                0.56, 0.58, 0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])[::-1] #[mag] array to store the B-V colours of the stars
#########################################################################################################################

# Now assume convection occurs for stars with T_eff < 8000 K (or 7870 K to be more precise based on the table)

def calculate_B_V(T_eff):
    '''
    This function takes the effective temperature of a star and returns the B-V colour of the star using the calibration from Böhm-Vitense (1970).
    The function assumes that the stars are of population type 1 and on the main sequence, and that convection occurs for stars with T_eff < 7870 K.

    Input: effective temperature [K]
    Output: corresponding B-V colour [mag]
    '''
    if T_eff < 7870:
        B_V = np.interp(T_eff, T_eff_convective_array, B_V_convective_array) 
    else:
        B_V = np.interp(T_eff, T_eff_radiative_array, B_V_radiative_array) 

    return B_V

def mass_to_B_V(mass):
    '''
    This function takes the mass of a star and returns the B-V colour of the star by first calculating the luminosity and radius of the star using the mass-luminosity and mass-radius relations, then calculating the effective temperature using the Stefan-Boltzmann law, and finally calculating the B-V colour using the calibration from Böhm-Vitense (1970).

    Input: stellar mass [kg]
    Output: corresponding B-V colour [mag]
    '''
    L = mass_to_luminosity(mass) #[W]
    R = mass_to_radius(mass) #[m]
    T_eff = calculate_T_eff(L, R) #[K]
    B_V = calculate_B_V(T_eff) #[mag]

    return B_V



