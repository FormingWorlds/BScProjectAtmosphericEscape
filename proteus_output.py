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
                dom_specie.append(f'{system}')
                current_inst.append(f'{inst}')
                current_mass.append(f'{m}')
                # T_inf.append()
                # massloss_flux_cgs.append()
                # massloss_si.append()
                # hom_alt.append()
                # exo_alt.append()
                # hom_pressure.append()
                # exo_pressure.append()
                # hom_T.append()
                # exo_T.append()
                
                
#(T_inf, flux_arr, flux_SI_arr, hom_height_arr, exo_height_arr, hom_pressure_arr, exo_pressure_arr, hom_temp_arr, exo_temp_arr)
                
export_dict = {
    'Dominant species' : dom_specie,
    'Instellation' : current_inst,
    'Planet mass' : current_mass
}
                
df = pd.DataFrame(export_dict)
df.to_csv('proteus_diffusion_limit_results.csv', index=False)