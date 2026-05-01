### Importing Proteus atmospheres that we want to examine ###
import pandas as pd

df_planet_properties = pd.read_csv("atmospheres/ProteusGoogleDrive/planet_bulk_properties_H2_atmospheres.csv", sep='\t')
#print(df_planet_properties)

case1 = df_planet_properties[df_planet_properties['Case'] == '1_F_earth'].iloc[0] #iloc needed to convert table to number. df[df..] mean we grab whole row where the case is true
case1000 = df_planet_properties[df_planet_properties['Case'] == '1000_F_earth'].iloc[0]

df_atmospheric_profile_case1 = pd.read_csv("atmospheres/ProteusGoogleDrive/H2_atmosphere_1_M_earth_1_F_earth.csv", sep='\t')
df_atmospheric_profile_case1000 = pd.read_csv("atmospheres/ProteusGoogleDrive/H2_atmosphere_1_M_earth_1000_F_earth.csv", sep='\t')
