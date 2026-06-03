import pandas as pd
import numpy as np
from scipy.constants import k, G, N_A


#define a class for easier&cleaner work
class Atmosphere:
    
    def __init__(self, temperature, pressure, height, Kzz, mmw, density, VMR, molar_masses, Natoms, Hfrac, diff_coeffs, planet_mass, planet_rad):
        """init function"""
        self.temperature = temperature  #[K]
        self.pressure    = pressure     #[Pa]
        self.height      = height       #[m]
        self.Kzz         = Kzz          #[cm^2/s]
        self.mmw         = mmw          #[kg/mol]
        self.density     = density      #[kg/m^3]
        self.VMR         = VMR          #[Volume mixing ratio]
        self.molar_masses = molar_masses #[kg/mol]
        self.Natoms      = Natoms       #amount of atoms in molecule
        self.Hfrac       = Hfrac        #frac of H atoms to total atoms in species
        self.diff_coeffs = diff_coeffs  #[cgs]
        self.planet_mass = planet_mass  #[kg]
        self.planet_rad  = planet_rad   #[m]
        
    
    def Kzz_extension(self):
        """extends Kzz up to 10^-13 bar assuming its a power law"""
        
        steps = 30
        
        K = self.Kzz
        p = self.pressure
        
        slope = (np.log10(K[-1]) - np.log10(K[-10]))/(np.log10(p[-1]) - np.log10(p[-10]))
        p_ext = np.logspace(np.log10(p[-1]), -8, steps)
        K_ext = ((p_ext/p[-1])**slope) * K[-1]
        
        K = np.concatenate((K, K_ext), axis=0)
        p = np.concatenate((p, p_ext), axis=0)
        
        return(K, p)
    
    def isothermal_extension(self):
        """extends T, VMR(const), MMW(const), rho and z isothermally up to 10^-13 bar"""
        
        steps = 30
        p = self.pressure
        p_ext = np.logspace(np.log10(p[-1]), -8, steps+1)[1:]  
        
        T = self.temperature
        VMR = self.VMR
        mmw = self.mmw
        rho = self.density
        z = self.height
        M = self.planet_mass
        R = self.planet_rad
        
        #extending the constant values: T, mmw, vmr
        T_ext = np.full(steps, T[-1])
        mmw_ext = np.full(steps, mmw[-1])
        T = np.concatenate((T, T_ext), axis=0)
        mmw = np.concatenate((mmw, mmw_ext), axis=0)
        
        for element in VMR:
            last_VMR = VMR[element][-1]
            VMR_ext = np.full(steps, last_VMR)
            VMR[element] = np.concatenate((VMR[element], VMR_ext), axis=0)
        
        #extending z (assuming constant g) and rho via HSE
        g_const = (G*M)/((R+z[-1])**2)
        H_arr = (k*T_ext)/((mmw_ext/N_A)*g_const)
        
        z_ext = z[-1] + (H_arr * np.log(p[-1]/p_ext))
        rho_ext = ((mmw_ext/N_A)*p_ext)/(k*T_ext)
        
        z = np.concatenate((z, z_ext), axis=0)
        rho = np.concatenate((rho, rho_ext), axis=0)
        
        return(T, VMR, mmw, rho, z)
        
        
        
def import_atmosphere(atmo_file, diff_file, system, planet_file, instellation):
    """makes an object of the class Atmosphere by extracting the needed information from a PROTEUS csv file (atmo_file), PROTEUS planet bulk parameter file (planet_file),
    and my diffusion coefficient file (diff_file). also needs 'system', which is a string of the main constituent, 'instellation' which is also a string from PROTEUS data. 
    makes all arrays aligned w increasing altitude"""
    
    df = pd.read_csv(atmo_file, delimiter='\t')
    df1 = pd.read_csv(diff_file)
    df2 = pd.read_csv(planet_file, delimiter='\t')
    
    atmo = Atmosphere(
        
        temperature = np.flip(df['Temperature [K]'].values),
        pressure    = np.flip(df['Pressure [Pa]'].values),
        height      = np.flip(df['Height [m]'].values),
        Kzz         = np.flip(df['Kzz [cm2/s]'].values),
        mmw         = np.flip((df['MMW [g/mol]'].values)) * 1e-3,
        density     = np.flip(df['Density [kg/m3]'].values),
        
        VMR = {
            "H2" : np.flip(df['H2 [VMR]'].values),
            "H2O": np.flip(df['H2O [VMR]'].values),
            "H"  : np.flip(df['H [VMR]'].values),
            "CH4": np.flip(df['CH4 [VMR]'].values),
            "NH3": np.flip(df['NH3 [VMR]'].values),
            "H2S": np.flip(df['H2S [VMR]'].values),
            "N2": np.flip(df['N2 [VMR]'].values),
            "CO2": np.flip(df['CO2 [VMR]'].values),
            "CO" : np.flip(df['CO [VMR]'].values),
            "O" : np.flip(df['O [VMR]'].values),
            "O2" : np.flip(df['O2 [VMR]'].values),
            "S" : np.flip(df['S [VMR]'].values),
            "SO2" : np.flip(df['SO2 [VMR]'].values),
            "N": np.flip(df['N [VMR]'].values)
        },
        
        molar_masses = {
            "H2" : 2.016*1e-3,
            "H2O": 18.015*1e-3,
            "H"  : 1.008*1e-3,
            "CH4": 16.043*1e-3,
            "NH3": 17.031*1e-3,
            "H2S": 34.08*1e-3,
            "N2": 28.012*1e-3,
            "CO2": 44.009*1e-3,
            "CO" : 28.010*1e-3,
            "O" : 15.999*1e-3,
            "O2" : 31.998*1e-3,
            "S" : 32.059*1e-3,
            "SO2" : 64.066*1e-3,
            "N" : 14.007*1e-3
        },
        
        Natoms = {
            "H2" : 2,
            "H2O": 3,
            "H"  : 1,
            "CH4": 5,
            "NH3": 4,
            "H2S": 3
        },
        
        Hfrac = {
            "H2" : 1,
            "H2O": 2/3,
            "H"  : 1,
            "CH4": 4/5,
            "NH3": 3/4,
            "H2S": 2/3
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
