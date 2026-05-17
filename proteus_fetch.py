import pandas as pd
import numpy as np

#define a class for easier&cleaner work
class Atmosphere:
    
    def __init__(self, temperature, pressure, height, Kzz, mmw, density, VMR, molar_masses, diff_coeffs, planet_mass, planet_rad):
        
        self.temperature = temperature  #[K]
        self.pressure    = pressure     #[Pa]
        self.height      = height       #[m]
        self.Kzz         = Kzz          #[cm^2/s]
        self.mmw         = mmw          #[kg/mol]
        self.density     = density      #[kg/m^3]
        self.VMR         = VMR          #[Volume mixing ratio]
        self.molar_masses = molar_masses #[kg/mol]
        self.diff_coeffs = diff_coeffs  #[cgs]
        self.planet_mass = planet_mass  #[kg]
        self.planet_rad  = planet_rad   #[m]
        
        
def import_atmosphere(atmo_file, diff_file, system, planet_file, instellation):
    """makes an object of the class Atmosphere by extracting the needed information from a PROTEUS csv file (atmo_file), PROTEUS planet bulk parameter file (planet_file),
    and my diffusion coefficient file (diff_file). also needs 'system', which is a string of the main constituent and 'instellation' which is also a string from PROTEUS data"""
    
    df = pd.read_csv(atmo_file, delimiter='\t')
    df1 = pd.read_csv(diff_file)
    df2 = pd.read_csv(planet_file, delimiter='\t')
    
    atmo = Atmosphere(
        
        temperature = df['Temperature [K]'].values,
        pressure    = df['Pressure [Pa]'].values,
        height      = df['Height [m]'].values,
        Kzz         = df['Kzz [cm2/s]'].values,
        mmw         = (df['MMW [g/mol]'].values) * 1e-3,
        density     = df['Density [kg/m3]'].values,
        
        VMR = {
            "H2" : df['H2 [VMR]'].values,
            "H2O": df['H2O [VMR]'].values,
            "H"  : df['H [VMR]'].values,
            "CH4": df['CH4 [VMR]'].values,
            "NH3": df['NH3 [VMR]'].values,
            "H2S": df['H2S [VMR]'].values
        },
        
        molar_masses = {
            "H2" : 2.016*1e-3,
            "H2O": 18.015*1e-3,
            "H"  : 1.008*1e-3,
            "CH4": 16.043*1e-3,
            "NH3": 17.031*1e-3,
            "H2S": 34.08*1e-3
        },
        
        diff_coeffs = {
            "H2": {
                "A": (df1[(df1["minor_const"] == "H2") & (df1["major_const"] == system)])["A"].values,
                "s": (df1[(df1["minor_const"] == "H2") & (df1["major_const"] == system)])["s"].values
            },
            "H2O": {
                "A": (df1[(df1["minor_const"] == "H2O") & (df1["major_const"] == system)])["A"].values,
                "s": (df1[(df1["minor_const"] == "H2O") & (df1["major_const"] == system)])["s"].values
            },
            "H": {
                "A": (df1[(df1["minor_const"] == "H") & (df1["major_const"] == system)])["A"].values,
                "s": (df1[(df1["minor_const"] == "H") & (df1["major_const"] == system)])["s"].values
            },            
            "CH4": {
                "A": (df1[(df1["minor_const"] == "CH4") & (df1["major_const"] == system)])["A"].values,
                "s": (df1[(df1["minor_const"] == "CH4") & (df1["major_const"] == system)])["s"].values
            },
            "NH3": {
                "A": (df1[(df1["minor_const"] == "NH3") & (df1["major_const"] == system)])["A"].values,
                "s": (df1[(df1["minor_const"] == "NH3") & (df1["major_const"] == system)])["s"].values
            },
            "H2S": {
                "A": (df1[(df1["minor_const"] == "H2S") & (df1["major_const"] == system)])["A"].values,
                "s": (df1[(df1["minor_const"] == "H2S") & (df1["major_const"] == system)])["s"].values
            }
        },
        
        planet_mass = df2[df2["Case"] == instellation]["M_planet [kg]"].values,
        planet_rad = df2[df2["Case"] == instellation]["R_int [m]"].values
    )
    
    return atmo 

