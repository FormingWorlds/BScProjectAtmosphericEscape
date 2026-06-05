### This file will give birth to different stellar types which will be used to investigate the effect of different stellar types on the escape regime and rate of an atmosphere. ###
import numpy as np
import scipy as sp

solar_mass = 1.989 * 10**30 #kg
solar_radius = 6.957 * 10**8 #m

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
        log_L = 2.028 * np.log10(mass) - 0.976
    elif 0.45 < mass_sm <= 0.72:
        log_L = 4.572 * np.log10(mass) - 0.102
    elif 0.72 < mass_sm <= 1.05:
        log_L = 5.743 * np.log10(mass) - 0.007
    elif 1.05 < mass_sm <= 2.4:
        log_L = 4.329 * np.log10(mass) + 0.010
    elif 2.4 < mass_sm <= 7:
        log_L = 3.967 * np.log10(mass) + 0.093
    elif 7 < mass_sm <= 31:
        log_L = 2.865 * np.log10(mass) + 1.105
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




        
