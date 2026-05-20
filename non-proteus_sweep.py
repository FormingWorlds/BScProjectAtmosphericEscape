import numpy as np
from matplotlib import pyplot as plt
plt.rc('text', usetex=True)
from improved_limit import improved_limiting_flux
from classical_limit import classical_limiting_flux

#this code is for performing parameter sweeps of the classical and improved limiting fluxes and plotting them


#set variables for the function:

temperatures = np.linspace(100, 400, 50)  #exobase temperatures [K]
T_0          = 100                        #temperature at the bottom of the model (p_0), [K]
#temperatures = np.linspace(500, 1500, 50)  #exobase temperatures [K]
#T_0          = 500                        #temperature at the bottom of the model (p_0), [K]

minor_const     = 'H'                                 #diffusing minor constituent, all caps string
major_const     = 'CO2'                                #major background constituent, all caps string
M_minor         = 1*1e-3                              #molar mass of minor constituent, [kg]
M_major         = 44*1e-3                             #molar mass of major constituent, [kg]

#minor_const     = ['H', 'H', 'H', 'H', 'He', 'He', 'He']                              #diffusing minor constituent, all caps string
#major_const     = ['air', 'CO2', 'N2', 'O2', 'N2', 'air', 'O2']                                #major background constituent, all caps string
#M_minor         = np.array([1*1e-3, 1*1e-3, 1*1e-3, 1*1e-3, 4*1e-3, 4*1e-3, 4*1e-3])                           #molar mass of minor constituent, [kg]
#M_major         = np.array([28.9647*1e-3, 44*1e-3, 28*1e-3, 32*1e-3, 28*1e-3, 28.9647*1e-3, 32*1e-3])                              #molar mass of major constituent, [kg]

mole_frac_homop = 1e-5                                #constant mole fraction at the homopause [unitless]
#mole_frac_homop = np.array([1e-6, 1e-5, 1e-4, 1e-3])  #constant mole fraction at the homopause [unitless]

#M_p = 7.99*5.97e24  #8 x Earth, [kg]
#_p = 1.875*6371000  #2 x Earth, [m]
M_p = 6.417e23  #mass of Mars, [kg]
R_p = 3389.5e3  #radius of Mars, [m]

z_0 = 80e3  #height above surface at bottom of the model, [m]

p_0     = 0.1   #pressure defining bottom of the model, [Pa]
#p_0     = 100   #pressure defining bottom of the model, [Pa]
p_inf   = 1e-6   #pressure defining roughly the exobase, [Pa]
p_steps = 10000  #atmosphere layers

K = 3e6                   #eddy diffusion coefficient [cm^2/s]
#K = np.logspace(4, 7, 4)  #eddy diffusion coefficient [cm^2/s]
alpha = -0.25         #from Yelle 



#outputting results:
colors = ["#D55E00", "#0072B2", "#CC79A7", "#332288", "#661100", "#009E73", "#F0E442", "#88CCEE"]
flux_list = []
for i in range(0, 7):    
    for temp in temperatures:
        results = improved_limiting_flux(temp, minor_const[i], major_const[i], M_minor[i], M_major[i], mole_frac_homop, K, M_p, R_p, z_0, p_0, p_inf, p_steps, T_0, alpha)
        
        if isinstance(results, str):
            print(results)
            break
        else:            
            flux_list.append(results[0])
            fluxes = np.asarray(flux_list, dtype='float64')
    
    fluxes = np.reshape(fluxes, (i+1, 50))      
    plt.plot(temperatures, fluxes[i,:], label=f'{minor_const[i]} in {major_const[i]}', c=colors[i])

plt.xlabel('Exobase temperature [K]')
plt.ylabel('Escape flux [cm$^{-2}$ $\cdot$ s$^{-1}$]')
plt.yscale('log')
plt.grid()
plt.legend()
plt.tight_layout()
plt.title('escape flux limit from superEarth-like-planet as function of exobase temp and composition (Xi=1e-5, K=3e6)')
plt.savefig('plots/supearth_comp_sweep.png', dpi=300, bbox_inches='tight')
plt.show()

# #outputting results:
# colors = ["#D55E00", "#0072B2", "#CC79A7", "#332288", "#661100", "#009E73", "#F0E442", "#88CCEE"]

# for i in range(0, 50):    
#     X_list = []
#     xi_list = []
#     for temp in temperatures:
#         results = improved_limiting_flux(temp, minor_const, major_const, M_minor, M_major, mole_frac_homop, K, M_p, R_p, z_0, p_0, p_inf, p_steps, T_0, alpha)
        
#         if isinstance(results, str):
#             print(results)
#             break
#         else:            
#             X_list.append(results[9])
#             xi_list.append(results[10])
#             X = np.asarray(X_list, dtype='float64')
#             xi = np.asarray(xi_list, dtype='float64')
        
#     plt.plot(X, xi, label=f'T={temp}')

# plt.xlabel('H mole fraction')
# plt.ylabel('height proxy')
# plt.yscale('log')
# plt.grid()
# plt.legend()
# plt.tight_layout()
# plt.title('variation of H mole frac with height as function of exobase temp')
# plt.savefig('plots/mars_yelle_2.png', dpi=300, bbox_inches='tight')
# plt.show()


#troubeshooting print output
print(f"homopause index: {results[2]}, exobase index: {results[3]}, Dmax: {results[4]:g}, H_0: {results[5]:g}, H_max: {results[6]:g}, mfp_0: {results[7]:g}, mfp_max: {results[8]:g}")