import numpy as np
from scipy.constants import G, k, N_A

# def total_M_atmos(height, rho_atm, R):
# 	"""Calculate the total mass of the atmosphere
# 	using the height and density profiles, 
# 	as well as the radius of the planet
# 	"""

# 	int_values = []
# 	for r, d in zip(height, rho_atm):
# 		int_val = d * (R+r)**2
# 		int_values.append(int_val)

# 	integral = np.trapz(int_values)

# 	M_atm = 4 * np.pi * integral

# 	return M_atm

def total_M_atmos(height, rho_atm, R):
	"""Calculate the total mass of the atmosphere
	using the height and density profiles, 
	as well as the radius of the planet
	"""
	rho_atm = np.array(rho_atm)
	height = np.array(height)
    
	M_atm = np.trapz(4 * np.pi * rho_atm * (R+height)**2, height)

	return M_atm




def scaleheight(T, R, z, M, mmw):
	"""Calculate scaleheight"""
	
	g = (G * M) / ((R + np.array(z))**2)

	h = (k * np.array(T)) / ((np.array(mmw)/N_A) * 10**(-3) * g)
	h = np.array(h)
	avg_h = np.mean(h)

	return h, avg_h
