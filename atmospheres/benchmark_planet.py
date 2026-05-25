# According to literature, the planet which benchmarks the rr-limited regime is hot jupiter planet HD 189733 b.
# Want to see if my code reproduces that the planet is within the regime.

import scipy as sp
import numpy as np
import sys
import os

# This adds the parent directory to the search path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from rr_escape import examine_atmosphere_for_rr_escape
from atmospheres.atmosphere_setting import Atmosphere

G = sp.constants.G #[m^3 kg^-1 s^-2]
m_p = sp.constants.m_p #[kg]
k_b = sp.constants.k #[J K^-1]

# Properties of  HD 189733 b

#Quantity    Value                          Unit         Description
T_wind       = 10**4                        #[K]         temperature of the escaping atmosphere
mu_wind      = 0.62                           #[dimless]   total mean molecular weight of the escaping atmosphere assumed to be entirely ionised
mu_plus_wind = 1.3                           #[dimless]   mean molecular weight of the ions in the escaping atmosphere assumed to be entirely ionised

M_p          =  1.123 * 1.899 * 10**(27)    #[kg]        planetary mass
F_xuv        = 18                           #[kg s^-3]    bolometric flux
nu_0         = 4.835981008048 * 10**15       #[Hz]       20ev in Hz as given in Murray-Clay et al 2009, they didnt specify the ionising freq for H2O specifically, it almost seemed like they used this value for everything
R_p          = 1.138 * 7.149 *10**7         #[m]         planetary radius
P_0          = 2000                         #[Pa]        pressure at the optical photosphere
T_eq         = np.full(5000, 1250)           #[K]         equilibrium temperature of the planet
mu_photo     = 2.5                            #[dimless]   mean molecular weight at the optical photosphere, assumed to be entirely molecular hydrogen

radii = np.linspace(R_p, 5*R_p, 5000) #[m] array of radii from the planetary radius to 10 times the planetary radius

pressures = P_0 * np.exp(G * M_p * mu_photo * m_p / (k_b * T_eq) * (1/radii - 1/R_p) ) #[Pa] pressure profile of the atmosphere based on the barometric formula, where P_0: pressure at the optical photosphere, G: gravitational constant, M_p: planetary mass, mu_photo: mean molecular weight at the optical photosphere, m_p: proton mass, k_b: Boltzmann constant, T_eq: equilibrium temperature of the planet, radii: array of radii from the planetary radius to 10 times the planetary radius

atm_HD189733b = Atmosphere(
    T_wind=T_wind, 
    mu_wind=mu_wind, 
    M_p=M_p, 
    F_xuv=F_xuv, 
    nu_0=nu_0, 
    mu_plus_wind=mu_plus_wind, 
    R_p=R_p, 
    pressures=pressures, 
    temperatures=T_eq, #using the constant T_eq as a starting array
    heights=radii - R_p,
    dominant_species='H2He'
)

examine_atmosphere_for_rr_escape(atm_HD189733b)
print('On wikipedia it says astronomers have observed the mass-loss rate for HD 189733 b to be around 10^9 - 10^11 g/s, which is 10^6 - 10^8 kg/s, so if my code is correct then the mass loss rate should be around that range and the planet should be within the RR-limited regime.')