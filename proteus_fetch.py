### Importing Proteus atmospheres that we want to examine ###
import pandas as pd
from atmospheres.atmosphere_setting import Atmosphere

df_planet_properties = pd.read_csv("atmospheres/ProteusGoogleDrive/planet_bulk_properties_H2_atmospheres.csv", sep='\t')
#print(df_planet_properties)

case1 = df_planet_properties[df_planet_properties['Case'] == '1_F_earth'].iloc[0] #iloc needed to convert table to number. df[df..] mean we grab whole row where the case is true
case1000 = df_planet_properties[df_planet_properties['Case'] == '1000_F_earth'].iloc[0]

df_atmospheric_profile_case1 = pd.read_csv("atmospheres/ProteusGoogleDrive/H2_atmosphere_1_M_earth_1_F_earth.csv", sep='\t')
df_atmospheric_profile_case1000 = pd.read_csv("atmospheres/ProteusGoogleDrive/H2_atmosphere_1_M_earth_1000_F_earth.csv", sep='\t')


atm_H2_case1 = Atmosphere(
    #bulk properties
    M_p=case1['M_planet [kg]'],
    R_p=case1['R_int [m]'],
    F_xuv=case1['F_ins [W/m2]'],

    #composition dependent properties
    T_wind=10**4, 
    nu_0 = 3.288467085473 * 10**15, 
    mu_wind=0.5, 
    mu_plus_wind=1,

    #atmospheric profiles
    temperatures=df_atmospheric_profile_case1['Temperature [K]'].values,
    pressures=df_atmospheric_profile_case1['Pressure [Pa]'].values,
    heights=df_atmospheric_profile_case1['Height [m]'].values
)

atm_H2_case1000 = Atmosphere(
    #bulk properties
    M_p=case1000['M_planet [kg]'],
    R_p=case1000['R_int [m]'],
    F_xuv=case1000['F_ins [W/m2]'],

    #composition dependent properties
    T_wind=10**4, 
    nu_0 = 3.288467085473 * 10**15, 
    mu_wind=0.5, 
    mu_plus_wind=1,

    #atmospheric profiles
    temperatures=df_atmospheric_profile_case1000['Temperature [K]'].values,
    pressures=df_atmospheric_profile_case1000['Pressure [Pa]'].values,
    heights=df_atmospheric_profile_case1000['Height [m]'].values
)