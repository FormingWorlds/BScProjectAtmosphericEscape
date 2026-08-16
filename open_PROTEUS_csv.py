import pandas as pd

def open_PROTEUS_csv(filename, sep = '\t'):
	"""Opens the PROTEUS atmosphere csv files using  Pandas.
	Returns density, height, temperature, mean molecular weight 
	of the atmosphere from surface going upwards
	"""

	file = pd.read_csv(filename, sep = sep)

	rho = file["Density [kg/m3]"].tolist()
	rho = rho[::-1]

	height = file["Height [m]"].tolist()
	height = height[::-1]

	T = file["Temperature [K]"].tolist()
	T = T[::-1]

	mmw = file["MMW [g/mol]"].tolist()
	mmw = mmw[::-1]
	
	return rho, height, T, mmw


def open_bulk_PROTEUS(filename, sep = '\t'):
	"""Opens the bulk properties files from PROTEUS
	Gives time, planet radius, planet mass, flux and MMW
	"""

	file = pd.read_csv(filename, sep = sep)

	mass = file["M_planet [kg]"].tolist()
	mass_1F = mass[1]
	mass_1000F = mass[0]

	radius = file["R_int [m]"].tolist()
	rad_1F = radius[1]
	rad_1000F = radius[0]

	return mass_1F, mass_1000F, rad_1F, rad_1000F



