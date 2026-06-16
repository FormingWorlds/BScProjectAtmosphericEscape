'''
Malina Ovesen

Radiation-recombination-limited escape mechanism.
'''

import numpy as np
import scipy as sp
from conversions import cm_to_m
from el_escape import el_escape_rate



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
    Calculates the radius to the sonic point of the atmosphere.

    Takes input parameters: planetary mass [kg], speed of sound for medium [m s^-1], and the radius to the XUV photosphere [m]

    All calculations done in SI units.
    '''
    ### calculates R_s based on input ###
    G = sp.constants.G #[m^3 kg^-1 s^-2] 
    R_s_calc = G * M_p / (2 * c_s**2) #[m] radius to the sonic point, where G: gravitational constant, M_p: planetary mass, c_s: sound speed
    # checks if R_s is smaller than R_base, if so, sets R_s = R_base and prints a message, otherwise keeps R_s = G * M_p / (2 * c_s**2) and prints a message
    if R_s_calc < R_base:
        #In this case it is not transonic wind. We set R_s equal to R_base
        R_s = R_base
        return R_s, False, R_s_calc #we return the sonic point, whether wind is transonic and what the calculated R_s was.
    else: #So if sonic point further out than R_base
        R_s = R_s_calc #we set R_s equal to the one we calculated
        return R_s, True, R_s_calc #returns the sonic point, whether wind is transonic, and what the calculated R_s was-


def calc_density_at_sonic_point(M_p, F_xuv, nu_0, R_s, c_s, R_base, mu_plus_wind, rr_coeff):
    '''
    Calculates the density at the sonic point.

    Takes input parameters: x [unit], y [unit], z [unit], ...

    All calculations done in SI units. So it will convert the rr_coeff from cm^3 s^-1 to m^3 s^-1 if it is given, and if not it will use the hydrogen formula for the recombination coefficient.
    '''

    h = sp.constants.h #[J s] 
    
    alpha_rec_B = rr_coeff * cm_to_m**3 #[m^3 s^-1] recombination coefficient for case B recombination for the dominant species in the escaping atmosphere, where rr_coeff: recombination coefficient for radiative case B recombination for the dominant species in the escaping atmosphere in cm^3 s^-1

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
 
    M_rr_rate = 4 * np.pi* rho_s * c_s * R_s**2 #[kg s^-1] radiation-recombination-limited escape rate REFERENCE: Lopez 2017 eq. 4, where rho_s: density at the sonic point, c_s: sound speed, R_s: radius to the sonic point

    return M_rr_rate 

def get_escape_diagnostics(atm):
    '''
    Calculates the radiation-recombination-limited escape rate and all the necessary parameters for the calculation.

    Takes input parameters: x [unit], y [unit], z [unit], ...

    All calculations done in SI units.
    '''
    #### RR escape rate ####
    c_s = calc_sound_speed(atm.T_wind, atm.mu_wind)
    R_s, is_transonic, R_s_calc = calc_sonic_point_radius(atm.M_p, c_s, atm.R_base)

    rho_s = calc_density_at_sonic_point(atm.M_p, atm.F_xuv, atm.nu_0, R_s, c_s, atm.R_base, atm.mu_plus_wind, atm.rr_coeff)
    escape_rate_rr = rr_escape_rate(rho_s, c_s, R_s)

    #### EL escape rate ####
    escape_rate_el = el_escape_rate(atm.epsilon_xuv, atm.F_xuv, atm.R_base, atm.M_p, K_tide=1)

    is_rr_limited = escape_rate_rr < escape_rate_el
    
    return {
        "R_base [m]": atm.R_base,
        "c_s [m/s]": c_s,
        "R_s [m]": R_s,
        "rho_s [kg/m^3]": rho_s,
        "escape_rate_rr [kg/s]": escape_rate_rr,
        "dominant_species": atm.dominant_species,
        "dominant_species_found_in_dict": atm.dominant_species_found_in_dict,
        "P_base_at_R_base [Pa]": atm.P_base_at_R_base,
        "is_transonic" : is_transonic,
        "R_s_calc [m]" : R_s_calc,
        "escape_rate_el [kg/s]": escape_rate_el,
        "is_rr_limited": is_rr_limited

    }

def examine_atmosphere_for_rr_escape(atm, P_base=10**(-4)):
    '''
    Examines whether the atmosphere is in the radiation-recombination-limited escape regime and prints the necessary diagnostics.

    Takes input parameters: x [unit], y [unit], z [unit], ...

    All calculations done in SI units.
    '''
    print("Examining the following atmosphere:")
    
    print(f"Dominant species at the base of the escaping atmosphere: {atm.dominant_species}")
    print(f"Dominant species used to set wind parameters: {atm.dominant_species_found_in_dict}")
    print(f"Bulk properties: M_p = {atm.M_p:.2g} kg, R_p = {atm.R_p:.2g} m, F_xuv = {atm.F_xuv:.2g} W/m^2")
    print(f"Wind properties: T_wind = {atm.T_wind:.2g} K, mu_wind = {atm.mu_wind:.2f}, nu_0 = {atm.nu_0:.2e} Hz, mu_plus_wind = {atm.mu_plus_wind:.2f}, rr_coeff = {atm.rr_coeff:.2g} cm^3 s^-1")
    print(f"Photosphere properties: P_0 = {atm.pressures[0]:.2g} Pa, T_eq = {atm.T[0]:.2g} K")
    print()

    ### compare with salz et al 2016 if can host hydrodynamic escape ###
    G = sp.constants.G #[m^3 kg^-1 s^-2]
    grav_pot_SI = - G * atm.M_p / atm.R_p #[m^2 s^-2] gravitational potential at the planetary radius, where G: gravitational constant, M_p: planetary mass, R_p: planetary radius
    grav_pot_cgs = grav_pot_SI * cm_to_m**(-2) #[cm^2 s^-2] gravitational potential at the planetary radius in cgs units, where grav_pot_SI: gravitational potential at the planetary radius in SI units
    grav_compare = np.log10(-grav_pot_cgs) #[log10(cm^2 s^-2)] logarithm of the gravitational potential at the planetary radius in cgs units, where grav_pot_cgs: gravitational potential at the planetary radius in cgs units
    salz_strong_grav_threshold = 13.6 #Larger than this value and gravity is too high for hydrodynamic escape. No EL and no wind for RR. REFERENCE Salz et. al. 2016
    salz_weak_grav_threshold = 13.11 
    if grav_compare < salz_weak_grav_threshold:
        print(f"Planet can host hydrodynamic escape through wind: {grav_compare} < {salz_weak_grav_threshold} (log10(grav pot [cm^2 s^-2])")
    elif grav_compare > salz_strong_grav_threshold:
        print(f"Planet cannot host hydrodynamic escape through wind: {grav_compare} > {salz_strong_grav_threshold} (log10(grav pot [cm^2 s^-2])")
    else:
        print(f"Planet is in intermediate regime of gravitationally binding atmosphere, wind strength declines rapidly")


    results = get_escape_diagnostics(atm)

    print()
    print("Escape diagnostics for the atmosphere:")
    print(f"Using values from this P_base: {results['P_base_at_R_base [Pa]']:.2g} Pa")
    print(f"R_base: {results['R_base [m]']:.2g} m")
    print(f"c_s: {results['c_s [m/s]']:.2g} m/s")
    print(f"R_s: {results['R_s [m]']:.2g} m")
    print(f"Is the wind transonic? {results['is_transonic']}")
    if results['is_transonic'] is False:
        print(f"The calculated sonic point radius was {results['R_s_calc [m]']:.2g} m")
    print(f"rho_s: {results['rho_s [kg/m^3]']:.2g} kg/m^3")
    print(f"RR Escape rate: {results['escape_rate_rr [kg/s]']:.2g} kg/s")
    print(f"EL Escape rate: {results['escape_rate_el [kg/s]']:.2g} kg/s")
    print(f"Is the escape rate RR limited? {results['is_rr_limited']}")
    
    

