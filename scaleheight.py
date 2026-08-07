import numpy as np

def scaleheight(z, rho):
	"""Calculate scaleheight"""
	lnrho = np.log(rho)
	b, a = np.polyfit(z, lnrho, 1)

	h = - 1.0 / b

	return h
