import numpy as np
from scipy.constants import G, k, N_A

# Constants used for scale height are imported by SciPy:
# gravitational constant G, Boltzmann constant k, Avogadro constant N_A


# Function to calculate the total mass of the atmosphere using density profile, by integrating using shells:

def total_M_atmos(height, rho_atm, R):
	"""Calculate the total mass of the atmosphere
	using the height and density profiles, 
	as well as the radius of the planet to integrate.
	"""
	rho_atm = np.array(rho_atm)
	height = np.array(height)
    
	M_atm = np.trapz(4 * np.pi * rho_atm * (R+height)**2, height)

	return M_atm


# Function to calculate the effective and average scale height of an atmosphere:

def scaleheight(T, R, z, M, mmw):
	"""Calculate scale height of the astmosphere
	   Using the temperature, height, and mean molecular weight profiles,
	   As well as the radius and mass of the planet.
	"""
	# Calculating the gravitational acceleration:
	g = (G * M) / ((R + np.array(z))**2)
	
	# Calculating the effective scale height and average scale height:
	h = (k * np.array(T)) / ((np.array(mmw)/N_A) * 10**(-3) * g)
	h = np.array(h)
	avg_h = np.mean(h)

	return h, avg_h
