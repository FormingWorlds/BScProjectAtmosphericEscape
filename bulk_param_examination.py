### Importing the escape functions ###
from rr_escape import examine_atmosphere_for_rr_escape, get_rr_escape_diagnostics
import scipy as sp
import numpy as np

### importing atmosphere class ###
from atmospheres.atmosphere_setting import Atmosphere

### Defining function to create atmosphere for loops ###
def make_simple_atmosphere(M_p, T_eq, mu_photo, F_xuv=None, F_ins=None,  P_0=2000, T_wind=10**4, dominant_species=None, resolution=5000, determine_radius=False, R_p=None, rr_coeff=None, P_base=0.0001,  nu_0=None, mu_wind=None, mu_plus_wind=None):
    '''
    Makes a simple isothermal atmosphere with the given input parameters.

    Takes input parameters: M_p [kg] - planetary mass, R_p [m] - planetary radius, F_ins [kg s^-3] - bolometric flux, nu_0 [Hz] - frequency of ionising radiation, mu_wind [dimless] - mean molecular weight of the escaping atmosphere, mu_plus_wind [dimless] - mean molecular weight of the ions in the escaping atmosphere, T_eq [K] - equilibrium temperature of the planet, mu_photo [dimless] - mean molecular weight at the optical photosphere, P_0 [Pa] - pressure at the optical photosphere.

    All calculations done in SI units.
    '''
    G = sp.constants.G #[m^3 kg^-1 s^-2]
    m_p = sp.constants.m_p #[kg]
    k_b = sp.constants.k #[J K^-1]

    if determine_radius:
        R_p = Atmosphere.determine_radius_from_MR_relation(M_p)
    elif determine_radius is False and R_p is not None:
        R_p = R_p
    else: 
        raise ValueError('You must either provide a planetary radius or set "determine_radius=True".') 
    radii = np.linspace(R_p, 10*R_p, resolution) #[m] array of radii from the planetary radius to 10 times the planetary radius

    pressures = P_0 * np.exp(G * M_p * mu_photo * m_p / (k_b * T_eq) * (1/radii - 1/R_p) ) #[Pa] pressure profile of the atmosphere based on the barometric formula, where P_0: pressure at the optical photosphere, G: gravitational constant, M_p: planetary mass, mu_photo: mean molecular weight at the optical photosphere, m_p: proton mass, k_b: Boltzmann constant, T_eq: equilibrium temperature of the planet, radii: array of radii from the planetary radius to 10 times the planetary radius

    return Atmosphere(
        T_wind=T_wind,
        mu_wind=mu_wind,
        M_p=M_p,
        F_xuv=F_xuv,
        F_ins=F_ins,
        P_base=P_base,
        nu_0=nu_0,
        mu_plus_wind=mu_plus_wind,
        R_p=R_p,
        pressures=pressures,
        temperatures=np.full(resolution, T_eq),
        heights=radii - R_p,
        dominant_species=dominant_species,
        rr_coeff=rr_coeff
    )


#want to loop over all the atmospheres with different planetary masses and get the escape diagnostics for each of them, and also find the break between rr regime not valid to valid
def sweep_parameter(atmospheres_varied_parameter, check_P_base=False, sensitivity=0.05, target_P_base=None):
    mass_loss_rates_varied_parameter = []
    regime_break_mass_varied_parameter = None
    rr_mask = []
    for i, atm in enumerate(atmospheres_varied_parameter):
        diagnostics = get_rr_escape_diagnostics(atm)
        mass_loss_rates_varied_parameter.append(diagnostics["escape_rate [kg/s]"])
        rr_mask.append(diagnostics["is_rr_limited"])
        if check_P_base:
            relative_diff = np.abs((atm.P_base_at_R_base - target_P_base[i])/target_P_base[i])
            if relative_diff > sensitivity:
                # Using a non-breaking warning print statement
                print(f"Boundary Mismatch at index {i} ({atm.dominant_species}): "
                      f"Used P_base ({atm.P_base_at_R_base:.4f} Pa) deviates from "
                      f"target P_base ({target_P_base} Pa) by {relative_diff:.4%} "
                      f"(Threshold: {sensitivity})")
        if diagnostics["R_s [m]"] > atm.R_base and regime_break_mass_varied_parameter is None:
            regime_break_mass_varied_parameter = i
    return mass_loss_rates_varied_parameter, regime_break_mass_varied_parameter, rr_mask
