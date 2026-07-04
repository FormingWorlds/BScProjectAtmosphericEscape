### This file will give birth to different stellar types which will be used to investigate the effect of different stellar types on the escape regime and rate of an atmosphere. ###
import numpy as np
import scipy as sp
from conversions import cm_to_m, erg_to_joule

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

    L = 10**log_L * solar_luminosity #[W]

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
#These temperatures correspond to K, G, F and A type stars. So mass-wise this means the input mass should be between around 0.45 solar masses and 2.1 solar masses based on the mass and T_eff of the stellar types.
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

def fetch_Jackson_parameters(B_V, original_mass=None):
    '''
    This function takes the B-V colour of a star and returns the parameters for the Jackson et. al. (2012) formula to calculate the X-ray-to-bolometric luminosity ratio.

    Input: B-V colour [mag], OPTIONAL: original stellar mass [kg]
    Output: log(L_X / L_bol) [Ø], log(tau_saturated) [yr], alpha (slope of the unsaturated regime) [Ø]
    '''
    if 0.290 <= B_V < 0.450:
        return -4.28, 7.87, 1.22
    
    elif 0.450 <= B_V < 0.565:
        return -4.24, 8.35, 1.24
    
    elif 0.565 <= B_V < 0.675:
        return -3.67, 7.84, 1.13
    
    elif 0.675 <= B_V < 0.790:
        return -3.71, 8.03, 1.28
    
    elif 0.790 <= B_V < 0.935:
        return -3.36, 7.90, 1.40
    
    elif 0.935 <= B_V < 1.275:
        return -3.35, 8.28, 1.09
    
    elif 1.275 <= B_V < 1.410:
        return -3.14, 8.21, 1.18
    else:
        mass_context = f" generated by a {original_mass/solar_mass:.4f} M_sun star" if original_mass is not None else ""
        raise ValueError(
            f"B-V value {B_V:.4f}{mass_context} is outside the bounds of Table 2 from Jackson et. al. (2012) (0.290 to 1.410)."
        )
    

def get_L_x_L_bol(B_V, age, original_mass=None):
    '''
    This function takes the B-V colour of a star and the age of the star and returns the X-ray luminosity of the star using the formula from Jackson et. al. (2012).

    Input: B-V colour [mag], age [yr], OPTIONAL: original stellar mass [kg]
    Output: X-ray luminosity [W]
    '''
    log_L_X_L_bol_saturated, log_tau_saturated, alpha = fetch_Jackson_parameters(B_V, original_mass)

    tau_saturated = 10**log_tau_saturated #[yr]
    L_X_L_bol_saturated = 10**log_L_X_L_bol_saturated

    if age < tau_saturated:
        L_X_L_bol = L_X_L_bol_saturated
    else:
        L_X_L_bol = L_X_L_bol_saturated * (age / tau_saturated)**(-alpha)

    return L_X_L_bol

def get_L_x(mass, age):
    '''
    This function takes the mass and age of a star and returns the X-ray luminosity of the star by first calculating the B-V colour using the mass-to-B-V function, then calculating the X-ray-to-bolometric luminosity ratio using the formula from Jackson et. al. (2012), and finally calculating the X-ray luminosity by multiplying the X-ray-to-bolometric luminosity ratio by the bolometric luminosity of the star.

    Input: stellar mass [kg], age [yr]
    Output: X-ray luminosity [W]
    '''
    B_V = mass_to_B_V(mass) #[mag]
    L_X_L_bol = get_L_x_L_bol(B_V, age) 
    L_bol = mass_to_luminosity(mass) #[W]

    L_X = L_X_L_bol * L_bol #[W]

    return L_X

def get_F(luminosity, semi_major_axis):
    '''
    This function calculates flux from luminosity and distance using the inverse square law.

    Input: luminosity of star [W], distance from star to planet [m]
    Output: X-ray flux at planet [W m^-2]
    '''
 
    F = luminosity / (4 * np.pi * semi_major_axis**2) #[W m^-2]

    return F

def get_L_EUV(L_x, semi_major_axis):
    '''
    This function takes the X-ray luminosity received and returns the EUV luminosity recevied using the relation from King et. al. (2018) and power law parameters alpha and gamma from boundary energy choice X_ray range 5-100 Å and EUV range 100-912 Å.

    Input: X-ray luminosity [W], distance from star to planet [m]
    Output: EUV luminosity [W]
    '''
    alpha_cgs = 650 #erg cm^-2 s^-1
    alpha_si = alpha_cgs * erg_to_joule * cm_to_m**(-2) #[W m^-2]
    gamma = -0.45

    EUV_over_X_fraction = alpha_si * (L_x / (4 * np.pi * semi_major_axis**2))**gamma #[Ø]
    L_EUV = EUV_over_X_fraction * L_x #[W]

    return L_EUV

def get_escape_fluxes(mass, age, semi_major_axis):
    '''
    This function takes the mass and age of a star and the distance from the star to the planet and returns the bolometric and XUV flux received by the planet.

    Input: stellar mass [kg], stellar age [yr], distance from star to planet [m]
    Output: bolometric flux at planet [W m^-2], XUV flux at planet [W m^-2]
    '''
    L_bol = mass_to_luminosity(mass) #[W]
    F_bol = get_F(L_bol, semi_major_axis) #[W m^-2]

    L_X = get_L_x(mass, age) #[W]
    L_EUV = get_L_EUV(L_X, semi_major_axis)
    L_XUV = L_X + L_EUV #[W]

    F_XUV = get_F(L_XUV, semi_major_axis) #[W m^-2]

    return F_bol, F_XUV