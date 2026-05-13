import pandas as pd
import numpy as np
from proteus_fetch import Atmosphere, import_atmosphere

atm_archetype = ["CO2", "H2", "H2O", "N2"]  
instellations = ["1_F_earth", "1000_F_earth"]
mass = ["1_M_earth", "10_M_earth"]


diff_file = 'diffusion_coefficients.csv'

for specie in atm_archetype:
    system = f'{specie}'
    
    for inst in instellations:
        instellation = f'{inst}'
        
        for m in mass:
            if((m=="1_M_earth") & (inst=="1_F_earth") & (specie=="H2O")):
                pass
            else:
                atmo_file = f'PROTEUS/atmos/{specie}_atmosphere_{m}_{inst}.csv'
                planet_file = f'PROTEUS/planet/planet_bulk_properties_{specie}_atmosphere_{m}.csv'
                
                atm = import_atmosphere(atmo_file, diff_file, system, planet_file, instellation)