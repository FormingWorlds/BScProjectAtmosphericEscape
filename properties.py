def total_M_atmos(height, rho_atm, R):
	"""Calculate the total mass of the atmosphere using the height and density profiles, as well as the radius of the planet"""	
	import numpy as np

	int_values = []
	for r, d in zip(height, rho_atm):
		int_val = d * (R+r)**2
		int_values.append(int_val)

	integral = np.trapz(int_values)

	M_atm = 4 * np.pi * integral

	return M_atm


