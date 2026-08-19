import numpy as np
from properties import scaleheight
from constants import k, N_A, G, rho_pl



def r_min(rho, rho_pl, h):
	"""Minimum radius of planetesimal impacts
	* rho = atmosphere density planet
	* rho_pl = density planetesimal impactor
	* h = scaleheight planet atmosphere
	"""
	r_min = (3*rho / rho_pl)**(1/3)*h
	return r_min


def r_cap(rho, rho_pl, h, R):
	""" Cap size of planetesimal impactor radius
	* rho = atm. density planet
	* rho_pl = density planetesimal impactor
	* h = scaleheight planet atm.
	* R = radius planet
	"""
	r_cap = (3*np.sqrt(2*np.pi)*rho / (4*rho_pl))**(1/3) * (h*R)**(1/2)
	return r_cap


def r_gi(h, R):
	"""Radius where giant impacts dominate
	* h = scaleheight atm
	* R = radius Earth
	"""
	r_gi = (2*h*(R**2))**(1/3)
	return r_gi




def cap_mass(rho_atm, h, R):
	"""Total cap mass for terrestial planets"""
	M_cap = 2 * np.pi * rho_atm * (h**2) * R
	return M_cap



def ejected_mass_planetesimal(M_atm, h, R, rho, rho_pl):
	"""Ejected mass by a planetesimal impactor (r < r_gi)
	Once r > r_cap, the cap mass is
	"""
	r_min_val = r_min(rho, rho_pl, h)
	r_cap_val = r_cap(rho, rho_pl, h, R)
	r_gi_val = r_gi(h, R)

	r_range = np.linspace(r_min_val, r_gi_val, 1000)
	M_ejected = []

	for r in r_range:

		if r < r_min_val:
			M_eject = 0.0
			print("No mass loss possible due to too small impactor size.")

		elif r >= r_min_val and r < r_cap_val:
			m_imp = rho_pl * (4/3) * np.pi * r**3
			M_eject = (r_min_val/(2*r)) * (1 - (r_min_val/r)**2) * m_imp

		elif r >= r_cap_val:
			M_eject = cap_mass(rho, h, R)

		M_ejected.append(M_eject)

	M_ejected = np.array(M_ejected)

	return M_ejected, r_range


# def ejected_mass_planetesimal(M_atm, h, R, rho, rho_pl):
# 	"""Ejected mass by a planetesimal impactor (r < r_gi)
# 	Once r > r_cap, the cap mass is
# 	"""
# 	r_min_val = r_min(rho, rho_pl, h)
# 	r_cap_val = r_cap(rho, rho_pl, h, R)
# 	r_gi_val = r_gi(h, R)

# 	r = np.linspace(r_min_val, r_gi_val, 1000)

# 	no_regime = (r < r_min_val)
# 	regime1 = (r >= r_min_val) and (r < r_cap_val)
# 	regime2 = (r >= r_cap_val) and (r < r_gi)
# 	regime3 = (r >= r_gi)

# 	if no_regime.any():
# 		print("Radius of impactor too small to eject atmosphere")

# 	M_eject[regime1] = (r_min_val/(2*r)) * (1 - (r_min_val/r)**2)
# 	M_eject[regime2] = cap_mass(rho, h, R)

# 	if regime3.any():
# 		print("Impactor big enough to eject all mass: global loss")


# 	return M_eject



