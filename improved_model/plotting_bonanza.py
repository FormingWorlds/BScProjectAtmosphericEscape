#species vs jeans 

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
plt.rc('text', usetex=True)


from cycler import cycler  # used to define color cycles
colorstyle = ["#D55E00", "#0072B2", "#CC79A7", "#009E73", "#332288", "#DDCC77", "#661100", "#AA4499", 
              "#44AA99", "#56B4E9", "#E69F00", "#F0E442", "#882255", "#999933", "#117733", "#88CCEE"]
plt.rcParams['axes.prop_cycle'] = cycler('color', colorstyle)

atm_archetype = ["CO2", "H2", "H2O", "N2"]  
instellations = ["1_F_earth", "1000_F_earth"]
mass = ["1_M_earth", "10_M_earth"]
T_infs = np.array([200, 300, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000])


diff_jeans = pd.read_csv('proteus_diff_vs_jeans.csv')

fig = plt.figure(figsize=(8,6))
gs = fig.add_gridspec(2, 2, hspace=0, wspace=0)
(ax1, ax2), (ax3, ax4) = gs.subplots(sharex='col', sharey='row')


colors = {
    'H2': "#D55E00",
    'H2O': "#CC79A7",
    'CO2': "#0072B2",
    'N2': "#009E73"
}
specie_markers = {
    '1_M_earth': '*',
    '10_M_earth': 'd'
}
flux_markers = {
    '1_F_earth': 40,
    '1000_F_earth': 120
}


for atm in atm_archetype:
    for inst in instellations:
        for m in mass:  
            for T in T_infs:
                
                row = diff_jeans[(diff_jeans['specie'] == atm) & (diff_jeans['instellation'] == inst) & (diff_jeans['planet_mass'] == m) & (diff_jeans['T_inf'] == T)]
                
                lim_diss = float(row['Lim_flux_disso[kg/s]'].values[0])      
                lim_H = float(row['Lim_flux_H[kg/s]'].values[0])    
                
                jeans_H = float(row['Jeans_flux_H[kg/s]'].values[0])    
                jeans_disso = float(row['Jeans_flux_disso[kg/s]'].values[0])    
                
                lim_diss_cgs = float(row['Lim_flux_disso[cgs]'].values[0])    
                lim_H_cgs = float(row['Lim_flux_H[cgs]'].values[0])    
                    
                if(atm=='CO2'):   
                    if((m=='1_M_earth') and (inst=='1_F_earth')):
                        ax1.scatter(
                            T, jeans_H,
                            color="#332288",
                            alpha=0.5,
                            marker = 'x',
                            s=20
                        )
                        ax1.scatter(
                            T, lim_H,
                            color=colors[atm],
                            alpha=0.5,
                            marker = specie_markers[m],
                            s=flux_markers[inst]
                        )
                    elif((m=='10_M_earth') and (inst=='1_F_earth')):
                        ax2.scatter(
                            T, jeans_H,
                            color="#332288",
                            alpha=0.5,
                            marker = 'x',
                            s=20
                        )
                        ax2.scatter(
                            T, lim_H,
                            color=colors[atm],
                            alpha=0.5,
                            marker = specie_markers[m],
                            s=flux_markers[inst]
                        )
                    elif((m=='1_M_earth') and (inst=='1000_F_earth')):
                        ax3.scatter(
                            T, jeans_H,
                            color="#332288",
                            alpha=0.5,
                            marker = 'x',
                            s=20
                        )
                        ax3.scatter(
                            T, lim_H,
                            color=colors[atm],
                            alpha=0.5,
                            marker = specie_markers[m],
                            s=flux_markers[inst]
                        )
                    elif((m=='10_M_earth') and (inst=='1000_F_earth')):
                        ax4.scatter(
                            T, jeans_H,
                            color="#332288",
                            alpha=0.5,
                            marker = 'x',
                            s=20
                        )
                        ax4.scatter(
                            T, lim_H,
                            color=colors[atm],
                            alpha=0.5,
                            marker = specie_markers[m],
                            s=flux_markers[inst]
                        )
                         
                                 
ax1.set_yscale('log')
ax2.set_yscale('log')
ax3.set_yscale('log')
ax4.set_yscale('log')


fig.supxlabel('Exobase extension limit [K]')
fig.supylabel('Diffusion limited flux [kg $\cdot$ s$^{-1}$]')


legend_handles2 = [
    Line2D([0], [0], color=colors['CO2'],
        label='$M_{\mathrm{p}} = 1M_{\mathrm{Earth}}$, $F = 1F_{\mathrm{Earth}}$', marker=specie_markers['1_M_earth'], linewidth=0, markersize=5, alpha=0.5),

    Line2D([0], [0], color=colors['CO2'],
        label='$M_{\mathrm{p}} = 1M_{\mathrm{Earth}}$, $F = 1000F_{\mathrm{Earth}}$', marker=specie_markers['1_M_earth'], linewidth=0, markersize=10, alpha=0.5),
    
    Line2D([0], [0], color=colors['CO2'],
        label='$M_{\mathrm{p}} = 10M_{\mathrm{Earth}}$, $F = 1F_{\mathrm{Earth}}$', marker=specie_markers['10_M_earth'], linewidth=0, markersize=5, alpha=0.5),
    
    Line2D([0], [0], color=colors['CO2'],
        label='$M_{\mathrm{p}} = 10M_{\mathrm{Earth}}$, $F = 1000F_{\mathrm{Earth}}$', marker=specie_markers['10_M_earth'], linewidth=0, markersize=10, alpha=0.5),
    
    Line2D([0], [0], color="#332288",
        label='Calculated Jeans escape flux', marker='x', linewidth=0, markersize=5, alpha=0.5)
]
legend2 = fig.legend(handles=legend_handles2, loc='outside upper center', ncols=3, bbox_to_anchor=(0.52, 0.995),
    borderaxespad=0.) 


plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.show() 