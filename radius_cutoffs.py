import numpy as np

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


