import numpy as np
from matplotlib import pyplot as plt
plt.rc('text', usetex=True)
from improved_limit import improved_limiting_flux
from classical_limit import classical_limiting_flux


#this code is for performing parameter sweeps of the classical and improved limiting fluxes
#variables of interest
temperatures    = np.linspace(100, 400, 50)  #[K]
minor_const     = 'H'
major_const     = 'N2'
M_minor         = 1*1e-3               #[kg]
M_major         = 28*1e-3              #[kg]
mole_frac_homop = 1e-5                 #[unitless]
K               = np.logspace(6, 7, 4)               #[cm^2/s]


flux_list = []
for i in range(0, 4):    
    for temp in temperatures:
        results = improved_limiting_flux(temp, minor_const, major_const, M_minor, M_major, mole_frac_homop, K[i])
        
        flux_list.append(results[0])
        fluxes = np.asarray(flux_list, dtype='float64')
    
    fluxes = np.reshape(fluxes, (i+1, 50))      
    plt.plot(temperatures, fluxes[i,:], label=f'K={K[i]:g}')

plt.xlabel('Exobase temperature [K]')
plt.ylabel('Escape flux [cm$^{-2}$ $\cdot$ s$^{-1}$]')
plt.yscale('log')
plt.grid()
plt.legend()
plt.tight_layout()
plt.title('H escape flux limit from N2 atmosphere (Mars-like planet) as function of exobase temp and eddy coeff')
plt.savefig('plots/H_in_N2_temp_eddy.png', dpi=300, bbox_inches='tight')
plt.show()