import pandas as pd
import numpy as np

df0 = pd.read_csv('analysis/proteus_diffusion_limit_results_full_disso.csv') #disso case
df1 = pd.read_csv('analysis/proteus_diffusion_limit_results_only_H.csv') #undisso case

dfJ0 = pd.read_csv('analysis/proteus_jeans_case_summary.csv') #undisso case
dfJ1 = pd.read_csv('analysis/proteus_jeans_case_summary_dissociation.csv') #disso case


atm_archetype = ["CO2", "H2", "H2O", "N2"]  
instellations = ["1_F_earth", "1000_F_earth"]
mass = ["1_M_earth", "10_M_earth"]
T_infs = np.array([200, 300, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000])

dom_specie = []
current_inst = []
current_mass = []
T_inf = []
lim_disso = []
lim_H = []
jeans_H = []
jeans_disso = []
lim_disso_cgs = []
lim_H_cgs = []

for atm in atm_archetype:
    for inst in instellations:
        for m in mass:
            for T in T_infs:
                dom_specie.append(atm)
                current_inst.append(inst)
                current_mass.append(m)
                T_inf.append(T)
                
                row0 = df0[(df0['Dominant_spec'] == atm) & (df0['Instellation'] == inst) & (df0['Planet_mass'] == m) & (df0['T_inf'] == T)]
                row1 = df1[(df1['Dominant_spec'] == atm) & (df1['Instellation'] == inst) & (df1['Planet_mass'] == m) & (df1['T_inf'] == T)]
                rowJ0 = dfJ0[(dfJ0['atmosphere_type'] == atm) & (dfJ0['flux_case'] == inst) & (dfJ0['mass_case'] == m) & (dfJ0['T_inf'] == T)]
                rowJ1 = dfJ1[(dfJ1['atmosphere_type'] == atm) & (dfJ1['flux_case'] == inst) & (dfJ1['mass_case'] == m) & (dfJ1['T_inf'] == T)]
                

                lim_disso.append(row0['Lim_flux[kg/s]'].values[0])                
                lim_H.append(row1['Lim_flux[kg/s]'].values[0])
                
                
                if(np.shape(rowJ0['weighted_mass_loss_kg_s'].values)[0] == 0):
                    jeans_H.append('0')
                else:
                    jeans_H.append(rowJ0['weighted_mass_loss_kg_s'].values[0])
                
                if(np.shape(rowJ1['weighted_mass_loss_kg_s'].values)[0] == 0):
                    jeans_disso.append('0')
                else:
                    jeans_disso.append(rowJ1['weighted_mass_loss_kg_s'].values[0])
                
                lim_disso_cgs.append(row0['Lim_flux[1/cm^2_s]'].values[0])
                lim_H_cgs.append(row1['Lim_flux[1/cm^2_s]'].values[0])


export_dict = {
    'specie' : dom_specie,
    'instellation' : current_inst,
    'planet_mass' : current_mass,
    'T_inf' : T_inf,
    'Lim_flux_disso[kg/s]' : lim_disso,
    'Lim_flux_H[kg/s]' : lim_H,
    'Jeans_flux_H[kg/s]' : jeans_H,
    'Jeans_flux_disso[kg/s]' : jeans_disso,
    'Lim_flux_disso[cgs]' : lim_disso_cgs,
    'Lim_flux_H[cgs]' : lim_H_cgs
}

df = pd.DataFrame(export_dict)
df.to_csv('proteus_diff_vs_jeans.csv', index=False)