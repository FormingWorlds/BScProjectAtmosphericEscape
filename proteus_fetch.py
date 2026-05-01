### Importing Proteus atmospheres that we want to examine ###
import pandas as pd

df_planet_properties = pd.read_csv("atmospheres/ProteusGoogleDrive/planet_bulk_properties_H2_atmospheres.csv", sep='\t')
#print(df_planet_properties)

case_1 = df_planet_properties[df_planet_properties['Case'] == '1_F_earth'].iloc[0]
case_1000 = df_planet_properties[df_planet_properties['Case'] == '1000_F_earth'].iloc[0]
print(case_1)
print(case_1000)