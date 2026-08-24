import numpy as np
from planetesimal_and_MLR import r_min, r_cap, r_gi

# Earth constants to check the thresholds of impactor radius for Earth-like planet:
R_Earth = 6378e3	# Earth radius [m]
h_Earth = 8.5e3		# scaleheight [m]
rho_atm_Earth = 1.293	# density Earth atmosphere at 0 height [kg/m^3]
rho_pl = 2000		# density planetesimal [kg/m^3]

r_min_Earth = r_min(rho_atm_Earth, h_Earth, rho_pl)
print(f"r_min = {r_min_Earth} m, which should be ~1 km")

r_cap_Earth = r_cap(rho_atm_Earth, h_Earth, R_Earth, rho_pl)
print(f"r_min = {r_cap_Earth} m, which should be ~25 km")

r_gi_Earth = r_gi(h_Earth, R_earth)
print(f"r_min = {r_gi_Earth} m, which should be ~1000 km")
