import numpy as np
from radius_cutoffs import *
from Earth_constants import *
from figure14 import r

def cap_mass(rho_atm, h):
	"""Total cap mass for terrestial planets"""
	M_cap = 2 * np.pi * rho_atm * (h**2) * R
	return M_cap


# m_max = ((np.pi * R) / (2*h))**(1/2) * M_cap

# m_min = 4*np.pi*rho_atm*(h**3)

# m_imp = np.linspace(m_min, m_max, 1000)


r_min = r_min(rho_atm, rho_pl, h)

def M_ejected(r_min, r):
	"""Mass ejected per radius"""
	M_ej_list = []
	for _ in r:
		M_ej = (r_min/(2*_)) * (1 - (r_min/_)**2)
		M_ej_list.append(M_ej)
	
	M_eject = np.array(M_ej_list)
	M_eject = np.abs(M_eject)
	return M_eject



