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


disso_fracs = {   #fraction of the initial molecules that are fully dissociated (the physics module takes into account the differing amount of atoms in different molecules)
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
            if(((m=="1_M_earth") & (inst=="1000_F_earth") & (specie=="H2"))):
                pass
            elif(((m=="10_M_earth") & (inst=="1000_F_earth") & (specie=="H2O"))):
                pass
            else:
                atmo_file = f'PROTEUS/atmos/{specie}_atmosphere_{m}_{inst}.csv'
                planet_file = f'PROTEUS/planet/planet_bulk_properties_{specie}_atmosphere_{m}.csv'
                
                atm = import_atmosphere(atmo_file, diff_file, system, planet_file, instellation)
                res = proteus_improved_limiting_flux(atm, system, disso_fracs)
                
                if res is None:
                    print(f'{system}, {m}, {instellation}: homopause not found - process terminated')
                else:
                    #print(f"flux from {system}, {m}, {instellation}: {res}")   
                    plt.plot(res[1], res[0], '*')
                    plt.title(f'{system}, {m}, {instellation}; all dissofracs = 1')  
                    plt.ylabel('Escape flux [cm$^{-2}$ $\cdot$ s$^{-1}$]')
                    plt.xlabel('Exobase temperature [K]')   
                    plt.grid()
                    plt.show()       
                    
                    # plt.plot(res[1],res[0]/100000, label="mfp, $[m]$")
                    # plt.plot(res[2],res[0]/100000, label="H, $[m]$")
                    # plt.title(f'{system}, {m}, {instellation}')
                    # plt.ylabel('p $[bar]$')
                    # plt.xlabel('Distance $[m]$')
                    # plt.xscale('log')
                    # plt.yscale('log')
                    # plt.gca().invert_yaxis()
                    # plt.grid()
                    # plt.legend()
                    # #plt.savefig('working_plots/proteus_part/DKmfpz_vs_p_Bates_extended_working.png', dpi=300, bbox_inches='tight')
                    # plt.show()
                
                #print(f"{system}, {m}, {instellation}: hom_id = {res[0]}")
                
                # a=0
                # for element in atm.VMR:
                #     if atm.VMR[element][-1] > 0.005:
                #         a += atm.VMR[element][-1]
                #         print(element)
                # print(f"for {system} {instellation} {m}, a={a}")
                
                
                # plt.plot(res[2]/100000, res[1], label='H')
                # plt.plot(res[2]/100000, res[0], label='mfp')
                # plt.plot(res[2]/100000, res[3], label='D')
                # plt.plot(res[2]/100000, res[4], label='K')
                # plt.xlabel('pressure $[bar]$')
                # plt.ylabel('$[cm^2/s]$ or $[m]$')
                # plt.yscale('log')
                # plt.xscale('log')
                # plt.title(f'{system}, {m}, {instellation}')
                # plt.gca().invert_xaxis()
                # plt.grid()
                # plt.legend()
                # #plt.savefig('working_plots/proteus_part/D_K_mfp_H_vs_z_issue.png', dpi=300, bbox_inches='tight')
                # plt.show()
                
                
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
                
                
                