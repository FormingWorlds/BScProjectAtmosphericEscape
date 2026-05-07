### Importing Proteus atmospheres that we want to examine ###
import pandas as pd
from atmospheres.atmosphere_setting import Atmosphere


def fetch_vmrs(df_atmospheric_profile):
    '''
    Fetches the volume mixing ratios (vmrs) from the atmospheric profile dataframe.

    Takes input parameters: df_atmospheric_profile [pandas dataframe] - dataframe containing the atmospheric profile, including vmrs.

    All calculations done in SI units.
    '''
    vmrs = {}
    for column in df_atmospheric_profile.columns:
        if 'VMR' in column:
            species = column.split()[0].strip() #assumes column name is in format "Species VMR [unit]"
            vmrs[species] = df_atmospheric_profile[column].values #array of vmr values for the species as a function of radius, where df_atmospheric_profile[column]: column of the dataframe containing the vmr values for the species
    return vmrs

#I want to fetch my 4 atmosphere cases for the 2 irradiation cases for the two masses that I am given
#I will do this in a loop to avoid repeating code
proteus_atmosphere_cases = {}


for comp in ['H2', 'H2O', 'CO2', 'N2']:
    for M in ['1_M_earth', '10_M_earth']:
        file_path1 = f"atmospheres/ProteusGoogleDrive/{comp}_atmospheres/planet_bulk_properties_{comp}_atmosphere_{M}.csv"
        
        try:
            df_planet_properties = pd.read_csv(file_path1, sep='\t')
            
            for case_name in ['1_F_earth', '1000_F_earth']:
                try:
                    case = df_planet_properties[df_planet_properties['Case'] == case_name].iloc[0]
                    profile_path = f"atmospheres/ProteusGoogleDrive/{comp}_atmospheres/{case_name}/{comp}_atmosphere_{M}_{case_name}.csv"
                    df_atmospheric_profile = pd.read_csv(profile_path, sep='\t')
                    
                    # Create Atmosphere object
                    proteus_atmosphere_cases[f"{comp}_{M}_{case_name}"] = Atmosphere(
                        M_p=case['M_planet [kg]'],
                        R_p=case['R_int [m]'],
                        F_xuv=case['F_xuv [W/m2]'], 
                        F_ins=case['F_ins [W/m2]'],
                        temperatures=df_atmospheric_profile['Temperature [K]'].values,
                        pressures=df_atmospheric_profile['Pressure [Pa]'].values,
                        heights=df_atmospheric_profile['Height [m]'].values,
                        vmrs=fetch_vmrs(df_atmospheric_profile)
                    )
                except FileNotFoundError:
                    print(f"Skipping profile: {comp}_{M}_{case_name} (File not found)")
                    continue  #makes it move onto the next step of the loop
                except IndexError:
                    print(f"Skipping case: {case_name} not found in bulk properties for {comp}_{M}")
                    continue

        except FileNotFoundError:
            print(f"Skipping group: Bulk properties for {comp}_{M} not found.")
            continue # Moves to the next M in the loop

