import numpy as np
from Earth_constants import *



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



def ejected_mass_planetesimal(r, M_atm, h, R, rho, rho_pl):
    """Ejected mass by a planetesimal impactor (r < r_gi)
       Once r > r_cap, the cap mass is 
    """
    r_min = r_min(rho, rho_pl, h)
    r_cap = r_cap(rho, rho_pl, h, R)
    
    if r < r_min:
        print("No mass loss possible due to escape.")
        
    elif r >= r_min and r < r_cap:
        M_eject = (r_min/(2*r)) * (1 - (r_min/r)**2)
    
    if r >= r_cap:
        M_eject = cap_mass(rho, h, R)
        
    return M_eject



