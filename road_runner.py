from rr_escape import *

### Importing the atmospheres we want to examine ###

from proteus_fetch import atm_H2_case1
from proteus_fetch import atm_H2_case1000
from atmospheres.simple_H2 import atm_H2
from atmospheres.simple_H2He_mix import atm_H2He_mix
from plots.plotter import plot_R_over_P, plot_R_over_T, plot_P_over_T
from scipy.constants import G


### Calculating the escape parameters ###

#I want to examine this atmosphere:
atm = atm_H2He_mix
print("Examining the following atmosphere:")
print(f"Bulk properties: M_p = {atm.M_p:.2g} kg, R_p = {atm.R_p:.2g} m, F_xuv = {atm.F_xuv:.2g} W/m^2")
print(f"Wind properties: T_wind = {atm.T_wind:.2g} K, mu_wind = {atm.mu_wind:.2f}, nu_0 = {atm.nu_0:.2e} Hz, mu_plus_wind = {atm.mu_plus_wind:.2f}")
print(f"Photosphere properties: P_0 = {atm.pressures[0]:.2g} Pa, T_eq = {atm.T[0]:.2g} K")
P_base = 10**(-4) #[Pa] pressure at the base of the escaping atmosphere, REFERENCE Lopez et. al. 2017
print()

### compare with salz et al 2016 if can host hydrodynamic escape ###
grav_pot_SI = - G * atm.M_p / atm.R_p #[m^2 s^-2] gravitational potential at the planetary radius, where G: gravitational constant, M_p: planetary mass, R_p: planetary radius
grav_pot_cgs = grav_pot_SI * cm_to_m**(-2) #[cm^2 s^-2] gravitational potential at the planetary radius in cgs units, where grav_pot_SI: gravitational potential at the planetary radius in SI units
grav_compare = np.log10(-grav_pot_cgs) #[log10(cm^2 s^-2)] logarithm of the gravitational potential at the planetary radius in cgs units, where grav_pot_cgs: gravitational potential at the planetary radius in cgs units
salz_strong_grav_threshold = 13.6 #Larger than this value and gravity is too high for hydrodynamic escape. No EL and no wind for RR. REFERENCE Salz et. al. 2016
salz_weak_grav_threshold = 13.11 #Smaller than this value and gravity is so weak the envelope will blowout. EL escape but no RR since no subsonic region. REFERENCE Salz et. al. 2016
print(f"Planet can host hydrodynamic (el) escape at least: {grav_compare < salz_strong_grav_threshold and grav_compare > salz_weak_grav_threshold} (log10(grav pot [cm^2 s^-2]) = {grav_compare:.2f}, where the strong gravity threshold is {salz_strong_grav_threshold} and the weak gravity threshold is {salz_weak_grav_threshold})")
print()

plot_R_over_P(atm)
plot_R_over_T(atm)
plot_P_over_T(atm)
results = get_rr_escape_diagnostics(atm, P_base)

print("Escape diagnostics for the atmosphere:")
print(f"R_base: {results['R_base [m]']:.2g} m")
print(f"c_s: {results['c_s [m/s]']:.2g} m/s")
print(f"R_s: {results['R_s [m]']:.2g} m")
print(f"rho_s: {results['rho_s [kg/m^3]']:.2g} kg/m^3")
print(f"Escape rate: {results['escape_rate [kg/s]']:.2g} kg/s")
if results['R_s [m]'] == results['R_base [m]']:
    print("NB: Escape is not radiation-recombination-limited.")