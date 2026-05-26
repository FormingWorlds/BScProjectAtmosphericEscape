import pandas as pd
import numpy as np
from proteus_fetch import Atmosphere, import_atmosphere
from proteus_physics import proteus_improved_limiting_flux
from matplotlib import pyplot as plt
plt.rc('text', usetex=True)


from cycler import cycler  # used to define color cycles
colorstyle = ["#D55E00", "#0072B2", "#CC79A7", "#009E73", "#332288", "#DDCC77", "#661100", "#AA4499", 
              "#44AA99", "#56B4E9", "#E69F00", "#F0E442", "#882255", "#999933", "#117733", "#88CCEE"]
plt.rcParams['axes.prop_cycle'] = cycler('color', colorstyle)



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

for inst in instellations:
    instellation = f'{inst}'
    
    for specie in atm_archetype:
        system = f'{specie}'
        
        for m in mass:
            if((m=="1_M_earth") & (inst=="1_F_earth") & (specie=="H2O")):
                pass
            elif((m=="1_M_earth") & (inst=="1000_F_earth") & (specie=="H2")):
                pass
            else:
                atmo_file = f'PROTEUS/atmos/{specie}_atmosphere_{m}_{inst}.csv'
                planet_file = f'PROTEUS/planet/planet_bulk_properties_{specie}_atmosphere_{m}.csv'
                
                atm = import_atmosphere(atmo_file, diff_file, system, planet_file, instellation)
                res = proteus_improved_limiting_flux(atm, system, disso_fracs)
                
                # a=0
                # for element in atm.VMR:
                #     if atm.VMR[element][-1] > 0.005:
                #         a += atm.VMR[element][-1]
                #         print(element)
                # print(f"for {system} {instellation} {m}, a={a}")
                
                
                print(np.shape(res[8]))
                
                plt.plot(res[4], res[7], label='H')
                plt.plot(res[4], res[8], label='mfp')
                plt.plot(res[4], res[9], label='D')
                plt.plot(res[4], res[5], label='K')
                plt.xlabel('Height $[km]$')
                plt.ylabel('$[cm^2/s]$')
                plt.yscale('log')
                plt.xscale('log')
                #plt.title(f'{system}, {m}, {instellation}')
                plt.grid()
                plt.legend()
                #plt.savefig('working_plots/proteus_part/D_and_K_vs_z_extended.png', dpi=300, bbox_inches='tight')
                plt.show()
                
                
                # plt.plot(res[0], res[4]/1000, label=f"{system}, {instellation}, {m}")
                # plt.xlabel('Temperature $[K]$')
                # plt.ylabel('Height $[km]$')
                # plt.yscale('log')
                #plt.xscale('log')
                #print(res)
                #print('aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa')
    #plt.gca().invert_yaxis()
    #plt.grid()
    #plt.legend()
    #plt.savefig('working_plots/proteus_part/T_vs_z_iso_extended.png', dpi=300, bbox_inches='tight')
    #plt.show()
                