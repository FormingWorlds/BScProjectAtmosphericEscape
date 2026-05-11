### Importing the escape functions ###
from rr_escape import *

### importing atmosphere class ###
from atmospheres.atmosphere_setting import Atmosphere

### Defining function to create atmosphere for loops ###
def make_simple_atmosphere(M_p, R_p, nu_0, mu_wind, mu_plus_wind, T_eq, mu_photo, F_xuv=None, F_ins=None,  P_0=2000, T_wind=10**4):
    '''
    Makes a simple isothermal atmosphere with the given input parameters.

    Takes input parameters: M_p [kg] - planetary mass, R_p [m] - planetary radius, F_ins [kg s^-3] - bolometric flux, nu_0 [Hz] - frequency of ionising radiation, mu_wind [dimless] - mean molecular weight of the escaping atmosphere, mu_plus_wind [dimless] - mean molecular weight of the ions in the escaping atmosphere, T_eq [K] - equilibrium temperature of the planet, mu_photo [dimless] - mean molecular weight at the optical photosphere, P_0 [Pa] - pressure at the optical photosphere.

    All calculations done in SI units.
    '''
    G = sp.constants.G #[m^3 kg^-1 s^-2]
    m_p = sp.constants.m_p #[kg]
    k_b = sp.constants.k #[J K^-1]

    radii = np.linspace(R_p, 40*R_p, 1000) #[m] array of radii from the planetary radius to 40 times the planetary radius

    pressures = P_0 * np.exp(G * M_p * mu_photo * m_p / (k_b * T_eq) * (1/radii - 1/R_p) ) #[Pa] pressure profile of the atmosphere based on the barometric formula, where P_0: pressure at the optical photosphere, G: gravitational constant, M_p: planetary mass, mu_photo: mean molecular weight at the optical photosphere, m_p: proton mass, k_b: Boltzmann constant, T_eq: equilibrium temperature of the planet, radii: array of radii from the planetary radius to 10 times the planetary radius

    return Atmosphere(
        T_wind=T_wind,
        mu_wind=mu_wind,
        M_p=M_p,
        F_xuv=F_xuv,
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
Mass_array = np.linspace(0.5, 30, 100, endpoint=True) * M_earth #[kg] array of planetary masses from 1 to 30 times the mass of the earth

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

#want to loop over all the atmospheres with different planetary masses and get the escape diagnostics for each of them, and also find the break between rr regime not valid to valid
mass_loss_rates_H2_varied_mass = []
regime_break_mass_H2_varied_mass = None
for i, atm in enumerate(atmospheres_H2_varied_mass):
    diagnostics = get_rr_escape_diagnostics(atm)
    mass_loss_rates_H2_varied_mass.append(diagnostics["escape_rate [kg/s]"])
    if diagnostics["R_s [m]"] > atm.R_base and regime_break_mass_H2_varied_mass is None:
        regime_break_mass_H2_varied_mass = i


plot_M_planet_over_M_dot(Mass_array, mass_loss_rates_H2_varied_mass, regime_break_index=regime_break_mass_H2_varied_mass)

### Want to examine how planetary radius affects mass loss rate ###
R_earth = 6.371 * 10**6            #[m]         radius of the earth
Radius_array = np.linspace(1, 5, 30, endpoint=True) * R_earth #[m] array of planetary radii from 1 to 5 times the radius of the earth

atmospheres_H2_varied_radius = np.array([make_simple_atmosphere(
    M_p=8.41 * 5.9722 * 10**24,
    R_p=R,
    F_ins=10**2.93 * erg_to_joule * cm_to_m**(-2) * 10**6,
    nu_0=3.288467085473 * 10**15,
    mu_wind=0.5,
    mu_plus_wind=1.0,
    T_eq=np.full(1000, 553),
    mu_photo=2,
    P_0=2000,
    T_wind=10**4
) for R in Radius_array])

#want to loop over all the atmospheres with different planetary radii and get the escape diagnostics for each of them, and also find the break between rr regime not valid to valid
mass_loss_rates_H2_varied_radius = []
regime_break_radius_H2_varied_radius = None
for i, atm in enumerate(atmospheres_H2_varied_radius):
    diagnostics = get_rr_escape_diagnostics(atm)
    mass_loss_rates_H2_varied_radius.append(diagnostics["escape_rate [kg/s]"])
    if diagnostics["R_s [m]"] > atm.R_base and regime_break_radius_H2_varied_radius is None:
        regime_break_radius_H2_varied_radius = i

plot_R_planet_over_M_dot(Radius_array, mass_loss_rates_H2_varied_radius, regime_break_index=regime_break_radius_H2_varied_radius)



### Want to examine how F_xuv affects planetary mass loss rate ###

### Want to examine how planetary radius affects mass loss rate ###
F_xuv_earth =   0.2196388835          #[W m^-2]  Standard Earth XUV radiation
F_xuv_array = np.linspace(0.5, 1000, 300, endpoint=True) * F_xuv_earth #[m] array of planetary radii from 1 to 5 times the radius of the earth

atmospheres_H2_varied_F_xuv = np.array([make_simple_atmosphere(
    M_p=8.41 * 5.9722 * 10**24,
    R_p=2.73 * 6.371 * 10**6,
    F_xuv=F,
    nu_0=3.288467085473 * 10**15,
    mu_wind=0.5,
    mu_plus_wind=1.0,
    T_eq=np.full(1000, 553),
    mu_photo=2,
    P_0=2000,
    T_wind=10**4
) for F in F_xuv_array])

#want to loop over all the atmospheres with different xuv fluxes received and get the escape diagnostics for each of them, and also find the break between rr regime not valid to valid
mass_loss_rates_H2_varied_F_xuv = []
regime_break_radius_H2_varied_F_xuv = None
for i, atm in enumerate(atmospheres_H2_varied_F_xuv):
    diagnostics = get_rr_escape_diagnostics(atm)
    mass_loss_rates_H2_varied_F_xuv.append(diagnostics["escape_rate [kg/s]"])
    if diagnostics["R_s [m]"] > atm.R_base and regime_break_radius_H2_varied_F_xuv is None:
        regime_break_radius_H2_varied_radius = i

plot_xuv_flux_over_M_dot(F_xuv_array, mass_loss_rates_H2_varied_F_xuv, regime_break_index=regime_break_radius_H2_varied_F_xuv, compare=True)

