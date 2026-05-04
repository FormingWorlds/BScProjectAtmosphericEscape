from rr_escape import *

### Importing the atmospheres we want to examine ###

from proteus_fetch import atm_H2_case1
from proteus_fetch import atm_H2_case1000
from atmospheres.simple_H2 import atm_H2
from plotter import plot_P_over_R, plot_T_over_R


### Calculating the escape parameters ###

#I want to examine this atmosphere:
atm = atm_H2_case1000
P_base = 10**(-4) #[Pa] pressure at the base of the escaping atmosphere, REFERENCE Lopez et. al. 2017

plot_P_over_R(atm) 
plot_T_over_R(atm)
results = get_rr_escape_diagnostics(atm, P_base)

print("Escape diagnostics for the atmosphere:")
print(f"R_base: {results['R_base [m]']:.2g} m")
print(f"c_s: {results['c_s [m/s]']:.2g} m/s")
print(f"R_s: {results['R_s [m]']:.2g} m")
print(f"rho_s: {results['rho_s [kg/m^3]']:.2g} kg/m^3")
print(f"Escape rate: {results['escape_rate [kg/s]']:.2g} kg/s")
if results['R_s [m]'] == results['R_base [m]']:
    print("NB: Escape is not radiation-recombination-limited.")