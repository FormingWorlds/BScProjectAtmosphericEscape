from rr_escape import *
from atmospheres.atmosphere_setting import Atmosphere
import atmospheres.simple_H2 as simple_H2
import atmospheres.simple_H2He_mix as simple_H2He_mix

### Importing the atmospheres we want to examine ###
atm_H2 = Atmosphere(
    T_wind=simple_H2.T_wind, 
    mu_wind=simple_H2.mu_wind, 
    M_p=simple_H2.M_p, 
    F_xuv=simple_H2.F_xuv, 
    nu_0=simple_H2.nu_0, 
    mu_plus_wind=simple_H2.mu_plus_wind, 
    R_p=simple_H2.R_p, 
    pressures=simple_H2.pressures, 
    temperatures=simple_H2.T_eq, #using the constant T_eq as a starting array
    heights=simple_H2.radii - simple_H2.R_p
)


### Calculating the escape parameters ###

#I want to examine this atmosphere:
atm = atm_H2
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