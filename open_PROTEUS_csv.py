def open_PROTEUS_csv(filename, sep = ","):
	"""Opens the PROTEUS atmosphere csv files using  Pandas.
	Returns density, height, temperature, mean molecular weight 
	of the atmosphere from surface going upwards
	"""
	import pandas as pd

	file = pd.read_csv(filename, sep = sep)

	rho = file["Density [kg/m3]"].tolist()
	rho = rho[::-1]

	height = file["Height [m]"].tolist()
	height = height[::-1]

	T = file["Temperature [K]"]
	T = T[::-1]

	mmw = file["MMW [g/mol]"]
	mmw = mmw[::-1]
	
	return rho, height, T, mmw

