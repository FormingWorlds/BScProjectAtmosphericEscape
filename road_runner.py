from rr_escape import *

### Importing the atmospheres we want to examine ###

from proteus_fetch import atm_H2_case1
from atmospheres.simple_H2 import atm_H2


### Calculating the escape parameters ###

#I want to examine this atmosphere:
atm = atm_H2_case1
P_base = 10**(-4) #[Pa] pressure at the base of the escaping atmosphere, REFERENCE Lopez et. al. 2017


R_base = find_R_base(P_base, atm.radii, atm.pressures)
print("R_base:", R_base, "m")

c_s = calc_sound_speed(atm.T_wind, atm.mu_wind)
print("c_s:", c_s, "m/s")

R_s = calc_sonic_point_radius(atm.M_p, c_s, R_base)
print("R_s:", R_s, "m")

rho_s = calc_density_at_sonic_point(atm.M_p, atm.F_xuv, atm.nu_0, atm.T_wind, R_s, c_s, R_base, atm.mu_plus_wind)
print("rho_s:", rho_s, "kg/m^3")

escape_rate = rr_escape_rate(rho_s, c_s, R_s)
print("escape_rate:", escape_rate, "kg/s")