import numpy as np
import matplotlib.pyplot as plt

R = 6378		# radius Earth [km]
M_planet = 5.98e24	# mass Earth [kg]
h = 8.85		# scaleheight Earth [km]
M_atm = 5.15e18		# Earth atmosphere mass [kg]
rho_0 = 1.293		# density at Earth surface [kg/m^3] (unused)
rho = 5.51e12		# density of Earth [kg/km^'3]

r_min = 1.0		# minimum impact size [km]
r_cap = 25.0		# cap size [km]

r_gi = (2 * h * (R**2))**(1/3) # where giant impacts dominate [km]
print(r_gi)		# giant impact radius, should be ~900 km

r = np.geomspace(1.000001, 3000, 1000) # set up x-axis scale (log)


def M_T(r, r_min, r_cap, R, h, M_atm, M_planet, rho):
	"""Calculate total impactor mass for different sized impactors"""
	r_gi = (2 * h * (R**2))**(1/3)

	if r <= r_cap:
		M_T = (2*r / r_min) * (1 - (r_min/r)**2)**(-1) * M_atm

	elif r_cap <= r <=r_gi:
		M_T = (4*np.pi/3) * rho * r**3 * (2*R/h)

	elif r > r_gi:
		M_T = 4 * M_planet

	MT_Mplanet = M_T / M_planet

	return MT_Mplanet


MT_Mp = [M_T(_, r_min, r_cap, R, h, M_atm, M_planet, rho) for _ in r]
MT_Mp = np.array(MT_Mp)

