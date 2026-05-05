### Importing the escape functions ###
from rr_escape import *

### importing atmosphere class ###
from atmospheres.atmosphere_setting import Atmosphere

### Defining function to create atmosphere for loops ###
def make_simple_atmosphere(M_p, R_p, F_ins, nu_0, mu_wind, mu_plus_wind, T_eq, mu_photo, P_0=2000, T_wind=10**4):
    '''
    Makes a simple isothermal atmosphere with the given input parameters.

    Takes input parameters: M_p [kg] - planetary mass, R_p [m] - planetary radius, F_ins [kg s^-3] - bolometric flux, nu_0 [Hz] - frequency of ionising radiation, mu_wind [dimless] - mean molecular weight of the escaping atmosphere, mu_plus_wind [dimless] - mean molecular weight of the ions in the escaping atmosphere, T_eq [K] - equilibrium temperature of the planet, mu_photo [dimless] - mean molecular weight at the optical photosphere, P_0 [Pa] - pressure at the optical photosphere.

    All calculations done in SI units.
    '''
    G = sp.constants.G #[m^3 kg^-1 s^-2]
    m_p = sp.constants.m_p #[kg]
    k_b = sp.constants.k #[J K^-1]

    radii = np.linspace(R_p, 5*R_p, 1000) #[m] array of radii from the planetary radius to 10 times the planetary radius

    pressures = P_0 * np.exp(G * M_p * mu_photo * m_p / (k_b * T_eq) * (1/radii - 1/R_p) ) #[Pa] pressure profile of the atmosphere based on the barometric formula, where P_0: pressure at the optical photosphere, G: gravitational constant, M_p: planetary mass, mu_photo: mean molecular weight at the optical photosphere, m_p: proton mass, k_b: Boltzmann constant, T_eq: equilibrium temperature of the planet, radii: array of radii from the planetary radius to 10 times the planetary radius

    return Atmosphere(
        T_wind=T_wind,
        mu_wind=mu_wind,
        M_p=M_p,
        F_ins=F_ins,
        nu_0=nu_0,
        mu_plus_wind=mu_plus_wind,
        R_p=R_p,
        pressures=pressures,
        temperatures=T_eq,
        heights=radii - R_p
    )


### Want to examine how planetary mass affects mass loss rate ###
M_earth = 5.9722 * 10**24               #[kg]        mass of the earth
Mass_array = np.linspace(1, 30, 30, endpoint=True) * M_earth #[kg] array of planetary masses from 1 to 30 times the mass of the earth

atmospheres_H2_varied_mass = np.array([make_simple_atmosphere(
    M_p=M,
    R_p=2.73 * 6.371 * 10**6,
    F_ins=10**2.93 * erg_to_joule * cm_to_m**(-2) * 10**6,
    nu_0=3.288467085473 * 10**15,
    mu_wind=0.5,
    mu_plus_wind=1.0,
    T_eq=np.full(1000, 553),
    mu_photo=2,
    P_0=2000,
    T_wind=10**4
) for M in Mass_array]) #[Atmosphere] array of Atmosphere objects with different planetary masses from 1 to 30 times the mass of the earth, where M_p: planetary mass

mass_loss_rates_H2_varied_mass = np.array([get_rr_escape_diagnostics(atm)["escape_rate [kg/s]"] for atm in atmospheres_H2_varied_mass]) #[kg/s] array of mass loss rates for the atmospheres with different planetary masses from 1 to 30 times the mass of the earth, where atm: Atmosphere object, get_rr_escape_diagnostics: function that calculates the escape diagnostics for a given atmosphere and returns a dictionary with the escape rate in kg/s


plot_M_planet_over_M_dot(Mass_array, mass_loss_rates_H2_varied_mass)

