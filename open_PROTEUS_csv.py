def open_PROTEUS_csv(filename, sep = ","):
	"""Opens the PROTEUS atmosphere csv files using  Pandas.
	Returns density and height of the atmosphere from surface going upwards
	"""
	import pandas as pd

	file = pd.read_csv(filename, sep = sep)

	rho = file["Density [kg/m3]"].tolist()
	rho = rho[::-1]
	height = file["Height [m]"].tolist()
	height = height[::-1]

	return rho, height

