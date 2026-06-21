import numpy as np
from matplotlib import pyplot as plt
plt.rc('text', usetex=True)
from improved_limit import improved_limiting_flux
from classical_limit import classical_limiting_flux

colors = ["#D55E00", "#0072B2", "#CC79A7", "#332288", "#009E73", "#661100", "#F0E442", "#88CCEE"]

#this code is for performing parameter sweeps of the classical and improved limiting fluxes and plotting them
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(2, 2)

ax1 = fig.add_subplot(gs[1, 0])
ax2 = fig.add_subplot(gs[1, 1])
ax3 = fig.add_subplot(gs[0, :])   #span both columns


#set variables for the function:

temperatures = np.linspace(100, 400, 50)  #exobase temperatures [K]
T_0          = 100                        #temperature at the bottom of the model (p_0), [K]

minor_const     = 'H'                                 #diffusing minor constituent, all caps string
major_const     = 'CO2'                                #major background constituent, all caps string
M_minor         = 1*1e-3                              #molar mass of minor constituent, [kg]
M_major         = 44*1e-3                             #molar mass of major constituent, [kg]

M_p = 6.417e23  #mass of Mars, [kg]
R_p = 3389.5e3  #radius of Mars, [m]

z_0 = 80e3  #height above surface at bottom of the model, [m]

p_0     = 0.1   #pressure defining bottom of the model, [Pa]
p_inf   = 1e-6   #pressure defining roughly the exobase, [Pa]
p_steps = 10000  #atmosphere layers

K = 3e6                   #eddy diffusion coefficient [cm^2/s]
alpha = -0.25         #from Yelle 

mole_frac_homop = np.array([1e-3, 1e-4, 1e-5, 1e-6])  #mole fraction at the homopause grid [unitless]


#outputting results:
flux_list = []
for i in range(0, 4):    
    for temp in temperatures:
        results = improved_limiting_flux(temp, minor_const, major_const, M_minor, M_major, mole_frac_homop[i], K, M_p, R_p, z_0, p_0, p_inf, p_steps, T_0, alpha)
        
        if isinstance(results, str):
            print(results)
            break
        else:            
            flux_list.append(results[0])
            fluxes = np.asarray(flux_list, dtype='float64')
    
    fluxes = np.reshape(fluxes, (i+1, 50))      
    ax1.plot(temperatures, fluxes[i,:], label=f'Mole frac.={mole_frac_homop[i]:.2g}', c=colors[i])

ax1.set_xlabel('Exobase temperature [K]', fontsize=13)
ax1.set_ylabel('Escape flux [cm$^{-2}$ $\cdot$ s$^{-1}$]', fontsize=13)
ax1.set_yscale('log')
ax1.legend(fontsize=11, loc='upper right')



K = np.logspace(4, 7, 4)  #eddy diffusion coefficient [cm^2/s]
mole_frac_homop = 1e-5    #constant mole fraction at the homopause [unitless]

#outputting results:
flux_list = []
for i in range(0, 4):    
    for temp in temperatures:
        results = improved_limiting_flux(temp, minor_const, major_const, M_minor, M_major, mole_frac_homop, K[i], M_p, R_p, z_0, p_0, p_inf, p_steps, T_0, alpha)
        
        if isinstance(results, str):
            print(results)
            break
        else:            
            flux_list.append(results[0])
            fluxes = np.asarray(flux_list, dtype='float64')
    
    fluxes = np.reshape(fluxes, (i+1, 50))      
    ax2.plot(temperatures, fluxes[i,:], label=f'K={K[i]:.2g} $cm^2/s$', c=colors[i])

ax2.set_xlabel('Exobase temperature [K]', fontsize=13)
ax2.set_ylabel('Escape flux [cm$^{-2}$ $\cdot$ s$^{-1}$]', fontsize=13)
ax2.set_yscale('log')
ax2.legend(fontsize=11, loc='lower left')



minor_const     = ['H', 'H', 'H', 'H2', 'He']                              #diffusing minor constituent, all caps string
major_const     = ['CO2', 'O2', 'N2', 'N2', 'N2']                                #major background constituent, all caps string
M_minor         = np.array([1*1e-3, 1*1e-3, 1*1e-3, 2*1e-3, 4*1e-3])                           #molar mass of minor constituent, [kg]
M_major         = np.array([44*1e-3, 32*1e-3, 28*1e-3, 28*1e-3, 28*1e-3])                              #molar mass of major constituent, [kg]

K = 3e6                   #eddy diffusion coefficient [cm^2/s]


flux_list = []
for i in range(0, 5):    
    for temp in temperatures:
        results = improved_limiting_flux(temp, minor_const[i], major_const[i], M_minor[i], M_major[i], mole_frac_homop, K, M_p, R_p, z_0, p_0, p_inf, p_steps, T_0, alpha)
        
        if isinstance(results, str):
            print(results)
            break
        else:            
            flux_list.append(results[0])
            fluxes = np.asarray(flux_list, dtype='float64')
    
    fluxes = np.reshape(fluxes, (i+1, 50))      
    ax3.plot(temperatures, fluxes[i,:], label=f'{minor_const[i]} in {major_const[i]}', c=colors[i])

ax3.set_xlabel('Exobase temperature [K]', fontsize=13)
ax3.set_ylabel('Escape flux [cm$^{-2}$ $\cdot$ s$^{-1}$]', fontsize=13)
ax3.set_yscale('log')
ax3.legend(fontsize=11)

ax1.tick_params(axis='both', which='major', labelsize=11)    
ax2.tick_params(axis='both', which='major', labelsize=11)    
ax3.tick_params(axis='both', which='major', labelsize=11)    

plt.tight_layout()
plt.show()