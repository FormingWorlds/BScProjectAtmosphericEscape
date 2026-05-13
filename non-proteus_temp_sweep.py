import numpy as np
from matplotlib import pyplot as plt
plt.rc('text', usetex=True)
from improved_limit import improved_limiting_flux
from classical_limit import classical_limiting_flux
from classical_limit_grad0 import classical_limiting_flux_grad0


#this code is for performing temperature sweeps and plots of the classical and improved limiting fluxes on test planets
#variables of interest
temperatures    = np.arange(100, 400)         #[K]
minor_const     = 'H'
major_const     = 'CO2'
M_minor         = 1*1e-3                      #[kg]
M_major         = 44*1e-3                     #[kg]
mole_frac_homop = 1e-5                        #[unitless]
K               = 3e6                         #[cm^2/s]


#reproducing Yelle fig1
#improved flux
flux_list_imp = []
T_list = []
flux_list_class = []
flux_list_class_grad0 = []

for temp in temperatures:
    results_imp = improved_limiting_flux(temp, minor_const, major_const, M_minor, M_major, mole_frac_homop, K)
    results_class = classical_limiting_flux(temp, minor_const, major_const, M_minor, M_major, mole_frac_homop, K)
    results_class_grad0 = classical_limiting_flux_grad0(temp, minor_const, major_const, M_minor, M_major, mole_frac_homop, K)
    
    flux_list_imp.append(results_imp[0])
    T_list.append(results_imp[1])
    flux_list_class.append(results_class)
    flux_list_class_grad0.append(results_class_grad0)

    fluxes_imp = np.asarray(flux_list_imp, dtype='float64')
    Temps = np.asarray(T_list, dtype='float64')
    fluxes_class = np.asarray(flux_list_class, dtype='float64')
    fluxes_class_grad0 = np.asarray(flux_list_class_grad0, dtype='float64')
    

plt.plot(T_list, fluxes_imp, label='improved limiting flux')
plt.plot(temperatures, fluxes_class, label='classical limiting flux')
plt.plot(temperatures, fluxes_class_grad0, label='classical limiting flux ignoring dT/dxi')
plt.yscale('log')
plt.title('escape fluxes as func of exospheric temp for Mars-like atmo')
plt.xlabel('Exobase temperature [K]')
plt.ylabel('Escape flux [cm$^{-2}$ $\cdot$ s$^{-1}$]')
plt.grid()
plt.legend()
plt.savefig('plots/flux_limit_comparison_mars.png', dpi=300, bbox_inches='tight')
plt.show()