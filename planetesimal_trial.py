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
r_min = 1
r_gi = (2*h*R)**(1/3)

m_small = 4*np.pi*rho_0*(h**3)
print(m_small)

m_large = rho_0*(np.pi*h*R)**(3/2)*2**(0.5)

#if m > m_small:

def M_T(r, r_min, r_cap, R, h, M_atm, rho, M_Earth):
"""describe"""
r_gi = (2*h*R)**(1/3)

if r<r_cap:
M_T = (2*r / r_min) * (1 - (r_min/r)**2)**(-1) * M_atm
elif r>= r_cap and r<r_gi:
M_T = 4*np.pi / 3 * rho * r**3 * (2*R / h)
elif r>= r_gi:
M_T = 4*M_Earth
