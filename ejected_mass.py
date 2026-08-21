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


def mass_loss_rate(r_range, q, M_atm, M_pl, h, R, rho, rho_pl):
	"""dM_atm/dt"""
	# Ejected mass by planetesimal impactors
	M_eject, r_range = ejected_mass_planetesimal(M_atm, h, R, rho, rho_pl)
	r_range = np.array(r_range)

	# mass of the impactor:
	m_pl = rho_pl * (4/3) * np.pi * (r_range)**3

	mass_ejec = (r)**(-q) * M_eject
	Mass_I = np.trapz(mass_int, r)

	N0_ejec = (r)**(-q) * m_pl
	N0_I = np.trapz(N0_ejec, r)

	# Mass loss rate in kg/s
	dM_dt = - M_pl * (Mass_I / N0_I)

	return dM_dt



