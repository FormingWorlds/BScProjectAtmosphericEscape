### This file will give birth to different stellar types which will be used to investigate the effect of different stellar types on the escape regime and rate of an atmosphere. ###
import numpy as np

def mass_to_luminosity(mass):
    '''
    This function takes in the mass of a star and returns the luminosity of the star in using the mass-luminosity relation from Eker et. al. 2018.
    Their luminosity relation consits of 6 piecewise equations which tie different ranges of stellar masses to a luminosity. As the paper states, the relation is only valid for stars between 0.179 and 31 solar masses, '
    therefore this function only covers that range as well.


    Input: stellar mass [kg]
    Output: stellar luminosity [W]
    '''

    #Relation uses solar masses so must convert before using the relation
    solar_mass = 1.989 * 10**30 #kg
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
        
