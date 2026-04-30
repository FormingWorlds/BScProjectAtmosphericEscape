import scipy as sp
import numpy as np

G = sp.constants.G #[m^3 kg^-1 s^-2]
m_p = sp.constants.m_p #[kg]
k_b = sp.constants.k #[J K^-1]

#A properties of a simple isothermal, single species (H_2) atmosphere

#Quantity    Value                          Unit         Description
T_wind       = 10**4                        #[K]         temperature of the escaping atmosphere
mu_wind      = 0.5                          #[dimless]   total mean molecular weight of the escaping atmosphere assumed to be entirely ionised
mu_plus_wind = 1.0                          #[dimless]   mean molecular weight of the ions in the escaping atmosphere assumed to be entirely ionised

#testcase values for the other necessary input parameters for the rr_escape_rate function using jupiter mass and radius, and a typical XUV flux for a close-in exoplanet
M_p          = 1.898 * 10**27               #[kg]        planetary mass
#R_base     = 7.1492 * 10**7                #[m]         radius of the base of the escaping atmosphere
F_xuv        = 1.0 * 10**9                  #[kg s^-3]    XUV flux
nu_0         = 3.288467085473 * 10**15      #[Hz]        frequency of Lyman-alpha radiation, which is the ionising radiation for HI
R_p          = 7.0 * 10**7                  #[m]         planetary radius
P_0          = 2000                         #[Pa]        pressure at the optical photosphere
T_eq         = 1500                         #[K]         equilibrium temperature of the planet
mu_photo     = 2                            #[dimless]   mean molecular weight at the optical photosphere, assumed to be entirely molecular hydrogen

radii = np.linspace(R_p, 5*R_p, 1000) #[m] array of radii from the planetary radius to 10 times the planetary radius

pressures = P_0 * np.exp(G * M_p * mu_photo * m_p / (k_b * T_eq) * (1/radii - 1/R_p) ) #[Pa] pressure profile of the atmosphere based on the barometric formula, where P_0: pressure at the optical photosphere, G: gravitational constant, M_p: planetary mass, mu_photo: mean molecular weight at the optical photosphere, m_p: proton mass, k_b: Boltzmann constant, T_eq: equilibrium temperature of the planet, radii: array of radii from the planetary radius to 10 times the planetary radius

