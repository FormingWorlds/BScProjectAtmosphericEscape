import pandas as pd
import numpy as np
from proteus_fetch import Atmosphere, import_atmosphere
from proteus_physics import proteus_improved_limiting_flux
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


fig, ax = plt.subplots(figsize=(8,6))

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
            res = proteus_improved_limiting_flux(atm, system, disso_fracs)
            
            if res is None:
                print(f'{system}, {m}, {instellation}: homopause not found - process terminated')
            else:
                #print(f"flux from {system}, {m}, {instellation}: {res}") 
                
                ax.scatter(
                    res[1], res[0],
                    color=colors[specie],
                    alpha=0.5,
                    marker = specie_markers[m],
                    s=80
                )
                                    
                #ax.set_xscale('log')
                ax.set_yscale('log')
                

                ax.set_xlabel('Exobase temperature [K]')
                ax.set_ylabel('Diffusion limited escape flux [cm$^{-2}$ $\cdot$ s$^{-1}$]')
                
                legend_handles1 = [
                    Line2D([0], [0], color=colors['H2O'],
                        label='H2O', marker='o', linewidth=0),
                    Line2D([0], [0], color=colors['H2'],
                        label='H2', marker='o', linewidth=0),
                    Line2D([0], [0], color=colors['CO2'],
                        label='CO2', marker='o', linewidth=0),
                    Line2D([0], [0], color=colors['N2'],
                        label='N2', marker='o', linewidth=0),
                ]
                legend1 = ax.legend(handles=legend_handles1, loc='upper left', title='Dominant species')
                ax.add_artist(legend1)
                
                legend_handles2 = [
                    Line2D([0], [0], color="#332288",
                        label='$M_{\mathrm{p}} = 1M_{\mathrm{Earth}}$', marker=specie_markers['1_M_earth'], linewidth=0),
                    
                    Line2D([0], [0], color="#332288",
                        label='$M_{\mathrm{p}} = 10M_{\mathrm{Earth}}$', marker=specie_markers['10_M_earth'], linewidth=0)
                ]
                legend2 = ax.legend(handles=legend_handles2, loc='upper center', title='Planet mass') 
                
                #ax.grid(True, which='both')


    #ax.invert_xaxis()
    plt.title(f'Limiting escape fluxes for all {inst} cases')
    plt.tight_layout()
    #plt.savefig('working_plots/proteus_part/all_DK_mfpH_mesopause_K_extended.png', dpi=300, bbox_inches='tight')
    plt.show() 
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
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
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                #     plt.plot(res[1],res[0]/100000, label="mfp, $[m]$")
                #     plt.plot(res[2],res[0]/100000, label="H, $[m]$")
                #     plt.title(f'{system}, {m}, {instellation}')
                #     plt.ylabel('p $[bar]$')
                #     plt.xlabel('Distance $[m]$')
                #     plt.xscale('log')
                #     plt.yscale('log')
                #     plt.gca().invert_yaxis()
                #     plt.grid()
                #     plt.legend()
                #     #plt.savefig('working_plots/proteus_part/DKmfpz_vs_p_Bates_extended_working.png', dpi=300, bbox_inches='tight')
                #     plt.show()
                
                # print(f"{system}, {m}, {instellation}: hom_id = {res[0]}")
           
           
           
           
           
           
           
           
           
                
    #             a=0
    #             for element in atm.VMR:
    #                 if atm.VMR[element][-1] > 0.005:
    #                     a += atm.VMR[element][-1]
    #                     print(element)
    #             print(f"for {system} {instellation} {m}, a={a}")
                
                
    #             plt.plot(res[2]/100000, res[1], label='H')
    #             plt.plot(res[2]/100000, res[0], label='mfp')
    #             plt.plot(res[2]/100000, res[3], label='D')
    #             plt.plot(res[2]/100000, res[4], label='K')
    #             plt.xlabel('pressure $[bar]$')
    #             plt.ylabel('$[cm^2/s]$ or $[m]$')
    #             plt.yscale('log')
    #             plt.xscale('log')
    #             plt.title(f'{system}, {m}, {instellation}')
    #             plt.gca().invert_xaxis()
    #             plt.grid()
    #             plt.legend()
    #             #plt.savefig('working_plots/proteus_part/D_K_mfp_H_vs_z_issue.png', dpi=300, bbox_inches='tight')
    #             plt.show()
                
                
    #             plt.plot(res[0], res[4]/1000, label=f"{system}, {instellation}, {m}")
    #             plt.xlabel('Temperature $[K]$')
    #             plt.ylabel('Height $[km]$')
    #             plt.yscale('log')
    #             plt.xscale('log')
    #             print(res)
    #             print('aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa')
    # plt.gca().invert_yaxis()
    # plt.grid()
    # plt.legend()
    # plt.savefig('working_plots/proteus_part/T_vs_z_iso_extended.png', dpi=300, bbox_inches='tight')
    # plt.show()
                
                
                
                
                
                
                
                
                      
                
                
                
                
                
                
                
                
                
                
         #K plotter       
                # plt.plot(res[0][:40]/100000, res[2][:40], label=f'{system}, {m}, {instellation}', alpha=0.7)
                # plt.xlabel('Pressure $[bar]$')
                # plt.ylabel('$[cm^2/s]$')
                # plt.yscale('log')
                # plt.xscale('log')
                # plt.gca().invert_xaxis()
                # plt.grid()
                # #plt.title(f'Eddy diff. coeff. profile of {system}, {m}, {instellation}')                
                # plt.tight_layout()
                # plt.savefig('working_plots/proteus_part/K_meso_extended.png', dpi=300, bbox_inches='tight')
                # plt.show()
                
                
                
            #if(((m=="1_M_earth") & (inst=="1000_F_earth") & (specie=="H2"))):
            #    pass
            #elif(((m=="10_M_earth") & (inst=="1000_F_earth") & (specie=="H2O"))):
            #    pass
            #else: