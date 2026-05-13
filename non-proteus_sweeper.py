import numpy as np
from matplotlib import pyplot as plt
plt.rc('text', usetex=True)
from improved_limit import improved_limiting_flux


#this code is for performing parameter sweeps of the classical and improved limiting fluxes

#variables of interest
temperatures    = np.arange(100, 400)  #[K]
minor_const     = 'H'
major_const     = 'CO2'
M_minor         = 1*1e-3               #[kg]
M_major         = 44*1e-3              #[kg]
mole_frac_homop = 1e-5                 #[unitless]
K               = 3e6                  #[cm^2/s]


#temperature sweep
flux_list = []
T_list = []

for temp in temperatures:
    results = improved_limiting_flux(temp, minor_const, major_const, M_minor, M_major, mole_frac_homop, K)
    flux_list.append(results[0])
    T_list.append(results[1])
    
fluxes = np.asarray(flux_list, dtype='float64')
temps = np.asarray(T_list, dtype='float64')

plt.plot(temps, fluxes)
plt.xlabel('Exobase temperature [K]')
plt.ylabel('Escape flux [cm$^{-2}$ $\cdot$ s$^{-1}$]')
plt.grid()
plt.savefig('plots/improved_lim.png')
plt.show()