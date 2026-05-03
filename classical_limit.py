import pandas as pd
import numpy as np
from scipy.constants import k, G, N_A
from scipy.integrate import quad


#Mars-like values from Jelle
#1% H in CO2 atmosphere, constant mixing ratio

T         = 100                                        #temperature profile, [K]
K         = 3e6                                        #eddy diffusion coefficient, [cm^2/sec]
alpha     = -0.25                                      #thermal diffusion factor
M_p       = 6.417e23                                   #mass of Mars, [kg]
Mars_rad  = 3389.5e3                                   #radius of Mars, [m]
M_CO2     = 44.0095*10e-3                              #molar mass, [kg/mol]
M_H       = 1*10e-3                                    #molar mass, [kg/mol]
mole_frac = 0.01                                       #mole fraction of H
m_a       = (mole_frac*M_H + (1-mole_frac)*M_CO2)/N_A  #mean molecular mass, [kg]
m_i       = M_H/N_A                                    #minor constituent mean molecular mass, [kg]
temp_grad = 0                                          #temperature gradient, [K/m]



#1. binary diffusion parameter, b

df = pd.read_csv("diffusion_coefficients.csv", comment="#")
#select a system here;
row = df[df["system"] == "H_in_CO2"].iloc[0]
A = row["A"]
s = row["s"]
b = A*(T**s) * 10e-2  #[1/m*sec]


#2. density profile

#we follow Yelle and a Martian example and set the bottom of the atmosphere, z_0 at 80 km, with p_0 = 0.1 Pa
#the top of the atmosphere is set at 500 km

p_0   = 0.1                                               #pressure at bottom of the atmosphere [Pa]
z_0   = 80e3                                              #bottom of atmosphere [m]
z_top = 500e3                                             #top of atmosphere, [m]
steps = 1000                                              #layers of the atmosphere for computation 
z     = np.linspace(z_0+Mars_rad, z_top+Mars_rad, steps)  #altitude grid, [m]


def density_scale_height(z, M_planet, T_z, temp_grad, ma_z):
    '''the integrable function of the density profile'''
    g_z = G*M_planet/(z**2)
    return((1/T_z)*(temp_grad) + (g_z*ma_z/k*T_z))


def density_profile(z, mu_z, T_z, p_0):
    '''calculates the density profile of the atmosphere in kg/m^3''' #currently for a constant T and mu profile
    rho_0 = (mu_z*p_0)/(k*T_z)
    
    integrals = []
    for i in range(1, np.shape(z)[0]):
        I, err = quad(density_scale_height, z[0], z[i], args=(M_p, T, temp_grad, m_a))
        integrals.append(I)
    
    ints = np.asarray(integrals, dtype='float64')
    return(rho_0 * np.e**(-ints))


print(density_profile(z, m_a, T, p_0))
    
