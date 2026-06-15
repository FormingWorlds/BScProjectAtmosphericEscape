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


diff_jeans = pd.read_csv('analysis/proteus_diff_vs_jeans.csv')

fig, ax = plt.subplots(figsize=(8,6))


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
                     
                if((atm=='N2')&(m=='1_M_earth')&(inst=='1000_F_earth')):
                    
                    ax.scatter(
                        T, lim_H_cgs,
                        color=colors[atm],
                        alpha=0.5,
                        marker = 'o',
                        s=80
                    )
                    ax.scatter(
                        T, lim_diss_cgs,
                        color=colors[atm],
                        alpha=0.5,
                        marker = '^',
                        s=80
                    )

                elif((atm=='H2O')&(m=='1_M_earth')&(inst=='1000_F_earth')):
                    ax.scatter(
                        T, lim_H_cgs,
                        color=colors[atm],
                        alpha=0.5,
                        marker = 'o',
                        s=80
                    )
                    ax.scatter(
                        T, lim_diss_cgs,
                        color=colors[atm],
                        alpha=0.5,
                        marker = '^',
                        s=80
                    )             
                                 
ax.set_yscale('log')

fig.supxlabel('Exobase extension limit [K]')
fig.supylabel('Diffusion limited flux [cm$^{-2}$ $\cdot$ s$^{-1}$]')

legend_handles1 = [
    Line2D([0], [0], color=colors['H2O'],
        label='H2O', marker='o', linewidth=0, alpha=0.5),
    Line2D([0], [0], color=colors['N2'],
        label='N2', marker='o', linewidth=0, alpha=0.5),
]
legend1 = fig.legend(handles=legend_handles1, loc='outside upper left', ncols=2, bbox_to_anchor=(0.102, 0.995),
    borderaxespad=0.)
fig.add_artist(legend1)

legend_handles2 = [
    Line2D([0], [0], color="#332288",
        label='Full dissociation', marker='^', linewidth=0, markersize=8, alpha=0.5),

    Line2D([0], [0], color="#332288",
        label='Only atomic H', marker='o', linewidth=0, markersize=8, alpha=0.5),

]
legend2 = fig.legend(handles=legend_handles2, loc='outside upper right', ncols=2, bbox_to_anchor=(0.98, 0.997),
    borderaxespad=0.) 


plt.tight_layout(rect=[0, 0, 1, 0.975])
plt.show() 