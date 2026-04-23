#A properties of a simple isothermal, single species (H_2) atmosphere

#Quantity    Value        Unit         Description
T_wind       = 10**4      #[K]         temperature of the escaping atmosphere
mu_wind      = 0.5        #[dimless]   total mean molecular weight of the escaping atmosphere assumed to be entirely ionised
mu_plus_wind = 1.0        #[dimless]   mean molecular weight of the ions in the escaping atmosphere assumed to be entirely ionised

#testcase values for the other necessary input parameters for the rr_escape_rate function using jupiter mass and radius, and a typical XUV flux for a close-in exoplanet
M_p = 1.898 * 10**27 #[kg]     planetary mass
R_base = 7.1492 * 10**7 #[m]   radius of the base of the escaping atmosphere
F_xuv = 1.0 * 10**9 #[W m^-2]  XUV flux
nu_0 = 2.46607 * 10**15 #[Hz]  frequency of Lyman-alpha radiation, which is the ionising radiation for HI