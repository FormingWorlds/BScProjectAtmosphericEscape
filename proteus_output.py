import pandas as pd
import numpy as np
from proteus_fetch import Atmosphere, import_atmosphere
from proteus_physics import proteus_improved_limiting_flux, disso_fracs

#this file runs the algorithm on all existing proteus files and returns a .csv file with relevant data pertaining to the runs
#all the tuning (like dissociation fractions) is done in proteus_physics.py


atm_archetype = ["CO2", "H2", "H2O", "N2"]  
instellations = ["1_F_earth", "1000_F_earth"]
mass = ["1_M_earth", "10_M_earth"]

diff_file = 'diffusion_coefficients.csv'

dom_specie = []
current_inst = []
current_mass = []
T_inf = []
massloss_flux_cgs = []
massloss_si = []
hom_alt = []
exo_alt = []
hom_pressure = []
exo_pressure = []
hom_T = []
exo_T = []

    
for specie in atm_archetype:
    system = f'{specie}'

    for inst in instellations:
        instellation = f'{inst}'
                
        for m in mass:
            atmo_file = f'PROTEUS/atmos/{specie}_atmosphere_{m}_{inst}.csv'
            planet_file = f'PROTEUS/planet/planet_bulk_properties_{specie}_atmosphere_{m}.csv'
            
            atm = import_atmosphere(atmo_file, diff_file, system, planet_file, instellation)
            res = proteus_improved_limiting_flux(atm, system, disso_fracs)
            
            if res is None:
                print(f'{system}, {m}, {instellation}: homopause not found - process terminated')
            else:  
                for i in range(0, np.shape(res[0])[0]):
                    dom_specie.append(f'{system}')
                    current_inst.append(f'{inst}')
                    current_mass.append(f'{m}')
                    T_inf.append(res[0][i])
                    massloss_flux_cgs.append(res[1][i])
                    massloss_si.append(res[2][0][i])
                    hom_alt.append(res[3][i])
                    exo_alt.append(res[4][i])
                    hom_pressure.append(res[5][i])
                    exo_pressure.append(res[6][i])
                    hom_T.append(res[7][i])
                    exo_T.append(res[8][i])
                
                
#export of results                    
export_dict = {
    'Dominant_spec' : dom_specie,
    'Instellation' : current_inst,
    'Planet_mass' : current_mass,
    'T_inf' : T_inf,
    'Lim_flux[1/cm^2_s]' : massloss_flux_cgs,
    'Lim_flux[kg/s]' : massloss_si,
    'Homopause_alt[km]' : hom_alt,
    'Exobase_alt[km]' : exo_alt,
    'Homopause_pressure[bar]' : hom_pressure,
    'Exobase_pressure[bar]' : exo_pressure,
    'Homopause_temp[K]' : hom_T,
    'Exobase_temp[K]' : exo_T
}
                
df = pd.DataFrame(export_dict)
df.to_csv('proteus_diffusion_limit_results_full_disso.csv', index=False)