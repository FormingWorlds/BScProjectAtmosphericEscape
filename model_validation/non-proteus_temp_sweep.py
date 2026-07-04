import numpy as np
from matplotlib import pyplot as plt
plt.rc('text', usetex=True)
from improved_limit import improved_limiting_flux
from classical_limit import classical_limiting_flux
from classical_limit_no_T_grad import classical_limiting_flux_grad0


#this code is for performing temperature sweeps and plots of the classical and improved limiting fluxes on test planets
#variables of interest
temperatures = np.arange(100, 400)  #[K]
T_0          = 100                  #temperature at the bottom of the model (p_0), [K]

z_0     = 80e3   #height above surface at bottom of the model, [m]
p_0     = 10000  #pressure defining bottom of the model, [Pa]
p_inf   = 1e-6   #pressure defining roughly the exobase, [Pa]
p_steps = 10000  #atmosphere layers

minor_const     = 'H'
major_const     = 'CO2'
M_minor         = 1*1e-3   #[kg]
M_major         = 44*1e-3  #[kg]
mole_frac_homop = 1e-5     #[unitless]
K               = 3e6      #[cm^2/s]
alpha           = -0.25    #from Yelle 

M_p = 6.417e23  #mass of Mars, [kg]
R_p = 3389.5e3  #radius of Mars, [m]


#reproducing Yelle fig1
#improved flux
flux_list_imp = []
T_list = []
flux_list_class = []
flux_list_class_grad0 = []
flux_list_yelle = []

for temp in temperatures:
    results_imp = improved_limiting_flux(temp, minor_const, major_const, M_minor, M_major, mole_frac_homop, K, M_p, R_p, z_0, p_0, p_inf, p_steps, T_0, alpha)
    results_class = classical_limiting_flux(temp, minor_const, major_const, M_minor, M_major, mole_frac_homop, K)
    results_class_grad0 = classical_limiting_flux_grad0(temp, minor_const, major_const, M_minor, M_major, mole_frac_homop, K)
    results_class_yelle = classical_limiting_flux(T_0, minor_const, major_const, M_minor, M_major, mole_frac_homop, K)
    
    flux_list_imp.append(results_imp[0])
    T_list.append(results_imp[1])
    flux_list_class.append(results_class)
    flux_list_class_grad0.append(results_class_grad0)
    flux_list_yelle.append(results_class_yelle)

    fluxes_imp = np.asarray(flux_list_imp, dtype='float64')
    Temps = np.asarray(T_list, dtype='float64')
    fluxes_class = np.asarray(flux_list_class, dtype='float64')
    fluxes_class_grad0 = np.asarray(flux_list_class_grad0, dtype='float64')
    fluxes_yelle = np.asarray(flux_list_yelle, dtype='float64')
    

plt.plot(temperatures, fluxes_yelle, label='Yelle-like classical limit', c="#009E73")
plt.plot(temperatures, fluxes_class, label='Classical limit', c="#0072B2")
plt.plot(T_list, fluxes_imp, label='Improved limit', c="#D55E00")

plt.yscale('log')
#plt.title('escape fluxes as func of exospheric temp for Mars-like atmo')
plt.xlabel('Exobase temperature [K]', fontsize=14)
plt.ylabel('Escape flux [cm$^{-2}$ $\cdot$ s$^{-1}$]', fontsize=14)
#plt.grid()
plt.legend(fontsize=12)
#plt.savefig('plots/flux_limit_comparison_mars.png', dpi=300, bbox_inches='tight')
plt.show()