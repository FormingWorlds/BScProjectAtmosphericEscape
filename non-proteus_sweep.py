import numpy as np
from matplotlib import pyplot as plt
plt.rc('text', usetex=True)
from improved_limit import improved_limiting_flux
from classical_limit import classical_limiting_flux

#this code is for performing parameter sweeps of the classical and improved limiting fluxes and plotting them


#set variables for the function:

temperatures = np.linspace(1500, 2500, 50)  #exobase temperatures [K]
T_0          = 1500                        #temperature at the bottom of the model (p_0), [K]

minor_const     = 'He'                                 #diffusing minor constituent, all caps string
major_const     = 'N2'                                #major background constituent, all caps string
M_minor         = 4*1e-3                              #molar mass of minor constituent, [kg]
M_major         = 28*1e-3                             #molar mass of major constituent, [kg]
#mole_frac_homop = 1e-5                                #constant mole fraction at the homopause [unitless]
mole_frac_homop = np.array([1e-6, 1e-4, 1e-3, 1e-3])  #constant mole fraction at the homopause [unitless]

M_p = 7.99*5.97e24  #mass of planet, [kg]
R_p = 1.875*6371000  #radius of planet, [m]

z_0 = 80e3  #height above surface at bottom of the model, [m]

p_0     = 10000    #pressure defining bottom of the model, [Pa]
p_inf   = 1e-6   #pressure defining roughly the exobase, [Pa]
p_steps = 10000  #atmosphere layers

K = 5e6                   #eddy diffusion coefficient [cm^2/s]
#K = np.logspace(6, 7, 4)  #eddy diffusion coefficient [cm^2/s]
alpha = -0.25         #from Yelle 



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
    plt.plot(temperatures, fluxes[i,:], label=f'molefrac={mole_frac_homop[i]:g}')

plt.xlabel('Exobase temperature [K]')
plt.ylabel('Escape flux [cm$^{-2}$ $\cdot$ s$^{-1}$]')
plt.yscale('log')
plt.grid()
plt.legend()
plt.tight_layout()
plt.title('He escape flux limit from N2 atmosphere (8 earth mass planet) as function of exobase temp and mole fraction')
plt.savefig('plots/supearth_He_in_N2_temp_molefrac.png', dpi=300, bbox_inches='tight')
plt.show()


#troubeshooting print output
print(f"homopause index: {results[2]}, exobase index: {results[3]}, Dmax: {results[4]:g}, H_0: {results[5]:g}, H_max: {results[6]:g}, mfp_0: {results[7]:g}, mfp_max: {results[8]:g}")
