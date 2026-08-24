import pandas as pd
import numpy as np
from properties import scaleheight, total_M_atmos

# Opens the PROTEUS files with profiles:
def open_PROTEUS_csv(filename, sep = '\t'):
	"""Opens the PROTEUS atmosphere csv files using  Pandas.
	Returns density, height, temperature, mean molecular weight 
	of the atmosphere from surface going upwards
	"""

	file = pd.read_csv(filename, sep = sep)

	return {
	    "rho": file["Density [kg/m3]"].iloc[::-1].to_numpy(),
	    "height": file["Height [m]"].iloc[::-1].to_numpy(),
	    "T": file["Temperature [K]"].iloc[::-1].to_numpy(),
	    "mmw": file["MMW [g/mol]"].iloc[::-1].to_numpy()
	}


# Opens the PROTEUS files with bulk data:
def open_bulk_PROTEUS(filename, sep = '\t'):
	"""Opens the bulk properties files from PROTEUS
	Gives time, planet radius, planet mass, flux and MMW
	"""

	file = pd.read_csv(filename, sep = sep)

	return {
	    "1_F": {
	        "mass": file["M_planet [kg]"].iloc[1],
	        "radius": file["R_int [m]"].iloc[1]
	    },
	    "1000_F": {
	        "mass": file["M_planet [kg]"].iloc[0],
	        "radius": file["R_int [m]"].iloc[0]
	    }
	}


# The data from the PROTEUS files:
# The data can be accessed through dictionary[element][mass][flux]

profiles = {}
bulk_data = {}

masses = ['1_M', '10_M']
fluxes = ['1_F', '1000_F']
elements = ['H2', 'H2O', 'N2', 'CO2']
colors = ['#006BA4', '#FF800E', '#ABABAB', '#595959']


# Creating the dictionary and adding PROTEUS values:
for e in elements:
    profiles[e] = {}
    bulk_data[e] = {}
    
    for m in masses:
        bulk_file = f"PROTEUS_data/Bulk_properties/planet_bulk_properties_{e}_atmospheres_{m}_earth.csv"
        bulk_data[e][m] = open_bulk_PROTEUS(bulk_file, sep="\t")
        
        profiles[e][m] = {}
        
        for f in fluxes:
            filename = f"PROTEUS_data/{e}_atmosphere_{m}_earth_{f}_earth.csv"
            profiles[e][m][f] = open_PROTEUS_csv(filename, sep="\t")
            
            
# Adding scale height and atmosphere mass to the dictionary
for e in elements:
    for m in masses:
        for f in fluxes:
            
            atm = profiles[e][m][f]
            bulk = bulk_data[e][m][f]
            
            atm_mass = total_M_atmos(atm["height"], atm["rho"], bulk["radius"])
            bulk_data[e][m][f]["atm_mass"] = atm_mass
            
            h_eff, h_avg = scaleheight(atm["T"], bulk["radius"], atm["height"], bulk["mass"], atm["mmw"])
            
            bulk_data[e][m][f]["h_avg"] = h_avg
            profiles[e][m][f]["h_eff"] = h_eff
            
                    