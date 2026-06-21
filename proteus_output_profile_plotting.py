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
                
                
                
                
                
                
    #             ax.scatter(
    #                 res[1], res[0],
    #                 color=colors[specie],
    #                 alpha=0.5,
    #                 marker = specie_markers[m],
    #                 s=80
    #             )
                                    
    #             #ax.set_xscale('log')
    #             ax.set_yscale('log')
                

    #             ax.set_xlabel('Exobase temperature [K]')
    #             ax.set_ylabel('Diffusion limited escape flux [cm$^{-2}$ $\cdot$ s$^{-1}$]')
                
    #             legend_handles1 = [
    #                 Line2D([0], [0], color=colors['H2O'],
    #                     label='H2O', marker='o', linewidth=0),
    #                 Line2D([0], [0], color=colors['H2'],
    #                     label='H2', marker='o', linewidth=0),
    #                 Line2D([0], [0], color=colors['CO2'],
    #                     label='CO2', marker='o', linewidth=0),
    #                 Line2D([0], [0], color=colors['N2'],
    #                     label='N2', marker='o', linewidth=0),
    #             ]
    #             legend1 = ax.legend(handles=legend_handles1, loc='upper left', title='Dominant species')
    #             ax.add_artist(legend1)
                
    #             legend_handles2 = [
    #                 Line2D([0], [0], color="#332288",
    #                     label='$M_{\mathrm{p}} = 1M_{\mathrm{Earth}}$', marker=specie_markers['1_M_earth'], linewidth=0),
                    
    #                 Line2D([0], [0], color="#332288",
    #                     label='$M_{\mathrm{p}} = 10M_{\mathrm{Earth}}$', marker=specie_markers['10_M_earth'], linewidth=0)
    #             ]
    #             legend2 = ax.legend(handles=legend_handles2, loc='upper center', title='Planet mass') 
                
    #             #ax.grid(True, which='both')


    # #ax.invert_xaxis()
    # plt.title(f'Limiting escape fluxes for all {inst} cases')
    # plt.tight_layout()
    # #plt.savefig('working_plots/proteus_part/all_DK_mfpH_mesopause_K_extended.png', dpi=300, bbox_inches='tight')
    # plt.show() 
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                #   ax1.plot(
                #     res[0]/100000, res[1],
                #     color=colors['D'],
                #     alpha=0.7,
                #     linestyle=linestyle[specie]
                # )
                
                # ax1.plot(
                #     res[0]/100000, res[2],
                #     color=colors['K'],
                #     alpha=0.7,
                #     linestyle=linestyle[specie]
                # )
                

                # # mfp and H on right axis
                # ax2.plot(
                #     res[0]/100000, res[3],
                #     color=colors['mfp'],
                #     alpha=0.7,
                #     linestyle=linestyle[specie]
                # )

                # ax2.plot(
                #     res[0]/100000, res[4],
                #     color=colors['H'],
                #     alpha=0.7,
                #     linestyle=linestyle[specie]
                # )     
                # ax1.axvline(x=1e-5/100000, color="#332288", linestyle="--", label="Isoth. extension begins", alpha=0.7)
                                    
                # ax1.set_xscale('log')
                # ax1.set_yscale('log')
                # ax2.set_yscale('log')
                

                # ax1.set_xlabel('Pressure [bar]')
                # ax1.set_ylabel('D, K $[cm^2 \ s{-1}]$')
                # ax2.set_ylabel('H, MFP $[cm]$')

                # legend_handles = [
                #     Line2D([0], [0], color=colors['D'],
                #         label='Binary diffusion coefficient'),

                #     Line2D([0], [0], color=colors['K'],
                #         label='Eddy diffusion coefficient'),

                #     Line2D([0], [0], color=colors['mfp'],
                #         label='H atom mean free path'),

                #     Line2D([0], [0], color=colors['H'],
                #         label='H atom scale height'),
                    
                #     Line2D([0], [0],
                #         color="#332288",
                #         linestyle="--",
                #         label='Isoth. extension begins')
                # ]
                # ax1.legend(handles=legend_handles, loc='lower right')
                # ax1.grid(True, which='both')


                # ax1.invert_xaxis()
                # plt.tight_layout()
                # #plt.savefig('working_plots/proteus_part/all_DK_mfpH_mesopause_K_extended.png', dpi=300, bbox_inches='tight')
                # plt.show()              
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                    
#                 plt.plot(res[1], res[0], '*', label=f'{system}, {m}, {instellation}')
#                 #plt.title(f'{system}, {m}, {instellation}; all dissofracs = 1')  
#                 plt.ylabel('Escape flux [cm$^{-2}$ $\cdot$ s$^{-1}$]')
#                 plt.xlabel('Exobase temperature [K]')   
# plt.yscale('log')
# plt.legend()
# plt.grid()
# plt.show()       
                    
                    
                    
                    
      
                
                
                

                
                
                
            #if(((m=="1_M_earth") & (inst=="1000_F_earth") & (specie=="H2"))):
            #    pass
            #elif(((m=="10_M_earth") & (inst=="1000_F_earth") & (specie=="H2O"))):
            #    pass
            #else:
            
            
            
            
            
            
#profile plotting