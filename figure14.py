import numpy as np
from radius_cutoffs import *
from Earth_constants import *


r_min = r_min(rho_atm, rho_pl, h)
r_cap = r_cap(rho_atm, rho_pl, h, R)
r_gi = r_gi(h,R)

r = np.geomspace((r_min + 0.0001), 3000000, 1000000) # set up x-axis scale (log)


def M_T(r, r_min, r_cap, r_gi, R, h, M_atm, M_planet, rho):
	"""Calculate total impactor mass for different sized impactors"""
	
	if r <= r_cap:
		M_T = (2*r / r_min) * (1 - (r_min/r)**2)**(-1) * M_atm

	elif r_cap <= r <=r_gi:
		M_T = (4*np.pi/3) * rho * r**3 * (2*R/h)

	elif r > r_gi:
		M_T = 4 * M_planet

	MT_Mplanet = M_T / M_planet

	return MT_Mplanet


MT_Mp = [M_T(_, r_min, r_cap, r_gi, R, h, M_atm, M, rho_pl) for _ in r]
MT_Mp = np.array(MT_Mp)
MT_Mp = np.abs(MT_Mp)

