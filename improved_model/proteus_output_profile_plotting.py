#p-T and D vs mfp

import pandas as pd
import numpy as np
from proteus_fetch import Atmosphere, import_atmosphere
from proteus_physics_for_profiles import proteus_profiles
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from scipy.signal import argrelmin
from proteus_physics import disso_fracs
plt.rc('text', usetex=True)


from cycler import cycler  # used to define color cycles
colorstyle = ["#D55E00", "#0072B2", "#CC79A7", "#009E73", "#332288", "#DDCC77", "#661100", "#AA4499", 
              "#44AA99", "#56B4E9", "#E69F00", "#F0E442", "#882255", "#999933", "#117733", "#88CCEE"]
plt.rcParams['axes.prop_cycle'] = cycler('color', colorstyle)



atm_archetype = ["CO2", "H2", "H2O", "N2"]  
instellations = ["1_F_earth", "1000_F_earth"]
mass = ["1_M_earth", "10_M_earth"]


diff_file = 'diffusion_coefficients.csv'


colors = {
    'H2': "#D55E00",
    'H2O': "#CC79A7",
    'CO2': "#0072B2",
    'N2': "#009E73"
}
specie_markers = {
    #'H2': 'P',
    '1_M_earth': '*',
    '10_M_earth': 'd',
    #'N2': 's'
}


for inst in instellations:
    instellation = f'{inst}'
    
    for specie in atm_archetype:
        system = f'{specie}'
        
        for m in mass:
            
            atmo_file = f'PROTEUS/atmos/{specie}_atmosphere_{m}_{inst}.csv'
            planet_file = f'PROTEUS/planet/planet_bulk_properties_{specie}_atmosphere_{m}.csv'
            
            atm = import_atmosphere(atmo_file, diff_file, system, planet_file, instellation)
            res = proteus_profiles(atm, system, disso_fracs)
            
            if res is None:
                print(f'{system}, {m}, {instellation}: homopause not found - process terminated')
            else:             
                #5 H 6 mfp  
                
                #print(np.shape(res[6]))
                
                homo = res[3][0]
                exo = int(res[4][0])
                
                fig, ax = plt.subplots(1, 2, figsize=(11,5))
                
                ax[1].plot(res[6][0][homo-1:exo+60], res[0][0][homo-1:exo+60]/100000, label=f'H mean free path; T_inf = {res[2][0]}K', alpha=0.9, c=colorstyle[0])
                ax[1].plot(res[5][0][homo-1:exo+60], res[0][0][homo-1:exo+60]/100000, label=f'H scale height; T_inf = {res[2][0]}K', alpha=0.9, c=colorstyle[0], linestyle='dashed')
                ax[1].axhline(y = res[0][0][int(res[4][0])]/100000, color =colorstyle[0], linestyle ="-.", alpha=0.3)
                
                ax[1].plot(res[6][9][homo-1:exo+60], res[0][9][homo-1:exo+60]/100000, label=f'H mean free path; T_inf = {res[2][9]}K', alpha=0.9, c=colorstyle[4])
                ax[1].plot(res[5][9][homo-1:exo+60], res[0][9][homo-1:exo+60]/100000, label=f'H scale height; T_inf = {res[2][9]}K', alpha=0.9, c=colorstyle[4], linestyle='dashed')
                ax[1].axhline(y = res[0][9][int(res[4][9])]/100000, color =colorstyle[4], linestyle ="-.", alpha=0.3)
                
                ax[1].set_yscale('log')
                ax[1].set_xscale('log')
                ax[1].yaxis.set_inverted(True)
                ax[1].set_xlabel('Scale height/mean free path [m]', fontsize=11)
                ax[1].set_ylabel('Pressure [bar]', fontsize=11)
                ax[1].legend(fontsize=11)
                
                 
                ax[0].plot(res[1][0], res[0][0]/100000, label=f'T_inf = {res[2][0]}K', alpha=0.6)
                ax[0].plot(res[1][0][int(res[4][0])], res[0][0][int(res[4][0])]/100000, 'x', c=colorstyle[0], markersize=6)
                ax[0].axhline(y = res[0][0][int(res[4][0])]/100000, color =colorstyle[0], linestyle ="-.", alpha=0.3)
                
                ax[0].plot(res[1][2], res[0][2]/100000, label=f'T_inf = {res[2][2]}K', alpha=0.6)
                ax[0].plot(res[1][2][int(res[4][2])], res[0][2][int(res[4][2])]/100000, 'x', c=colorstyle[1], markersize=6)
                
                ax[0].plot(res[1][4], res[0][4]/100000, label=f'T_inf = {res[2][4]}K', alpha=0.6)
                ax[0].plot(res[1][4][int(res[4][4])], res[0][4][int(res[4][4])]/100000, 'x', c=colorstyle[2], markersize=6)
                
                ax[0].plot(res[1][7], res[0][7]/100000, label=f'T_inf = {res[2][7]}K', alpha=0.6)
                ax[0].plot(res[1][7][int(res[4][7])], res[0][7][int(res[4][7])]/100000, 'x', c=colorstyle[3], markersize=6)
                
                ax[0].plot(res[1][9], res[0][9]/100000, label=f'T_inf = {res[2][9]}K', alpha=0.6)
                ax[0].plot(res[1][9][int(res[4][9])], res[0][9][int(res[4][9])]/100000, 'x', c=colorstyle[4], markersize=6, label='Exobase')
                ax[0].axhline(y = res[0][9][int(res[4][9])]/100000, color =colorstyle[4], linestyle ="-.", alpha=0.3)

                
                ax[0].axhline(y = res[0][0][res[3]]/100000, color ="#882255", linestyle ="--", label='Homopause', alpha=0.8)
                
                #fig.suptitle(f'p-T profiles of {system}, {m}, {instellation}')
                ax[0].set_xlabel('Temperature [K]', fontsize=11)
                ax[0].set_ylabel('Pressure [bar]', fontsize=11)
                ax[0].set_yscale('log')
                #plt.xscale('log')
                ax[0].yaxis.set_inverted(True)
                ax[0].legend(fontsize=11)
                ax[0].tick_params(axis='both', which='major', labelsize=11)    
                ax[1].tick_params(axis='both', which='major', labelsize=11)    
                
                plt.tight_layout()
                #plt.savefig(f'{specie}_{m}_{inst}_p-T_mfpH.png')
                #plt.savefig(f'{specie}_{m}_{inst}_p-T_mfpH.pdf')
                plt.show()
                plt.close() 