### Importing the escape functions ###
from rr_escape import *

### importing atmosphere class ###
from atmospheres.atmosphere_setting import Atmosphere

### Defining function to create atmosphere for loops ###
def make_simple_atmosphere(M_p, R_p, nu_0, mu_wind, mu_plus_wind, T_eq, mu_photo, F_xuv=None, F_ins=None,  P_0=2000, T_wind=10**4, dominant_species=None, resolution=5000):
    '''
    Makes a simple isothermal atmosphere with the given input parameters.

    Takes input parameters: M_p [kg] - planetary mass, R_p [m] - planetary radius, F_ins [kg s^-3] - bolometric flux, nu_0 [Hz] - frequency of ionising radiation, mu_wind [dimless] - mean molecular weight of the escaping atmosphere, mu_plus_wind [dimless] - mean molecular weight of the ions in the escaping atmosphere, T_eq [K] - equilibrium temperature of the planet, mu_photo [dimless] - mean molecular weight at the optical photosphere, P_0 [Pa] - pressure at the optical photosphere.

    All calculations done in SI units.
    '''
    G = sp.constants.G #[m^3 kg^-1 s^-2]
    m_p = sp.constants.m_p #[kg]
    k_b = sp.constants.k #[J K^-1]

    radii = np.linspace(R_p, 10*R_p, resolution) #[m] array of radii from the planetary radius to 10 times the planetary radius

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
        temperatures=np.full(resolution, T_eq),
        heights=radii - R_p,
        dominant_species=dominant_species
    )


#want to loop over all the atmospheres with different planetary masses and get the escape diagnostics for each of them, and also find the break between rr regime not valid to valid
def sweep_masses(atmospheres_varied_mass):
    mass_loss_rates_varied_mass = []
    regime_break_mass_varied_mass = None
    for i, atm in enumerate(atmospheres_varied_mass):
        diagnostics = get_rr_escape_diagnostics(atm)
        mass_loss_rates_varied_mass.append(diagnostics["escape_rate [kg/s]"])
        if diagnostics["R_s [m]"] > atm.R_base and regime_break_mass_varied_mass is None:
            regime_break_mass_varied_mass = i
    
    return mass_loss_rates_varied_mass, regime_break_mass_varied_mass


