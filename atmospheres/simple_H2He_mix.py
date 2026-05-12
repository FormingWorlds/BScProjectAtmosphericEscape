import scipy as sp
import numpy as np
from atmospheres.atmosphere_setting import Atmosphere
from conversions import *

G = sp.constants.G #[m^3 kg^-1 s^-2]
m_p = sp.constants.m_p #[kg]
k_b = sp.constants.k #[J K^-1]

#A properties of a simple isothermal, single species (H_2) atmosphere

#Quantity    Value                          Unit         Description
T_wind       = 10**4                        #[K]         temperature of the escaping atmosphere
mu_wind      = 0.62                          #[dimless]   total mean molecular weight of the escaping atmosphere assumed to be entirely ionised
mu_plus_wind = 1.3                          #[dimless]   mean molecular weight of the ions in the escaping atmosphere assumed to be entirely ionised

#testcase values for the other necessary input parameters for the rr_escape_rate function using jupiter mass and radius, and a typical XUV flux for a close-in exoplanet
M_p          = 8.41 * 5.9722 * 10**24               #[kg]        planetary mass
#R_base     = 7.1492 * 10**7                #[m]         radius of the base of the escaping atmosphere
F_xuv        = 10**2.93 * erg_to_joule * cm_to_m**(-2)                #[kg s^-3]    XUV flux. 100 times the solar flux of earth
nu_0         = 4.835981008048 * 10**15      #[Hz]        20eV. Maybe this value because of helium in mix. For ionising the atoms.
R_p          = 2.73 * 6.371 * 10**6            #[m]         planetary radius. 2 times earth radius
P_0          = 2000                         #[Pa]        pressure at the optical photosphere
T_eq         = np.full(5000, 553)                         #[K]         equilibrium temperature of the planet
mu_photo     = 2.5                            #[dimless]   mean molecular weight at the optical photosphere, assumed to be entirely molecular hydrogen

radii = np.linspace(R_p, 5*R_p, 5000) #[m] array of radii from the planetary radius to 10 times the planetary radius

pressures = P_0 * np.exp(G * M_p * mu_photo * m_p / (k_b * T_eq) * (1/radii - 1/R_p) ) #[Pa] pressure profile of the atmosphere based on the barometric formula, where P_0: pressure at the optical photosphere, G: gravitational constant, M_p: planetary mass, mu_photo: mean molecular weight at the optical photosphere, m_p: proton mass, k_b: Boltzmann constant, T_eq: equilibrium temperature of the planet, radii: array of radii from the planetary radius to 10 times the planetary radius

atm_H2He_mix = Atmosphere(
    T_wind=T_wind, 
    mu_wind=mu_wind, 
    M_p=M_p, 
    dominant_species='H2He',
    F_xuv=F_xuv, 
    nu_0=nu_0, 
    mu_plus_wind=mu_plus_wind, 
    R_p=R_p, 
    pressures=pressures, 
    temperatures=T_eq, #using the constant T_eq as a starting array
    heights=radii - R_p,
    )
