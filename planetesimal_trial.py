import numpy as np
import matplotlib.pyplot as plt

r_min = 1	# minimum impact (km)
rho_0 = 1.293	# density at earth (kg/m^3)
h =  8.5	# scale height earth atmosphere (km)
R = 12756	# Radius earth (km)
M_Earth = 5.97e24	# mass earth (kg)
M_atm = 5.15

r = np.linspace(1,1000, 1000)

r_cap = 25
r_min

m_small = 4*np.pi*rho_0*(h**3)
print(m_small)

m_large = rho_0*(np.pi*h*R)**(3/2)*2**(0.5)

#if m > m_small:


