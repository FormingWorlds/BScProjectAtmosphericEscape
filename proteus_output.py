import pandas as pd
import numpy as np
from proteus_fetch import Atmosphere, import_atmosphere
from proteus_physics import proteus_improved_limiting_flux
from matplotlib import pyplot as plt
plt.rc('text', usetex=True)


atm_archetype = ["CO2", "H2", "H2O", "N2"]  
instellations = ["1_F_earth", "1000_F_earth"]
mass = ["1_M_earth", "10_M_earth"]


diff_file = 'diffusion_coefficients.csv'


disso_fracs = {
    "H2" : 1,
    "H2O": 1,
    "H"  : 1,
    "CH4": 1,
    "NH3": 1,
    "H2S": 1
}


for specie in atm_archetype:
    system = f'{specie}'
    
    for inst in instellations:
        instellation = f'{inst}'
        
        for m in mass:
            if((m=="1_M_earth") & (inst=="1_F_earth") & (specie=="H2O")):
                pass
            else:
                atmo_file = f'PROTEUS/atmos/{specie}_atmosphere_{m}_{inst}.csv'
                planet_file = f'PROTEUS/planet/planet_bulk_properties_{specie}_atmosphere_{m}.csv'
                
                atm = import_atmosphere(atmo_file, diff_file, system, planet_file, instellation)
                
                res = proteus_improved_limiting_flux(atm, system, disso_fracs)
                
                plt.plot(res[0], res[1], label='K')
                plt.plot(res[0], res[2], label='D')
                plt.xlabel('height[m]')
                plt.ylabel('coeffs [cm^2/s]')
                #plt.yscale('log')
                plt.title(f'{system}, {m}, {instellation}')
                plt.grid()
                plt.legend()
                plt.show()