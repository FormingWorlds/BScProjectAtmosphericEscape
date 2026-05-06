import pandas as pd
import numpy as np
from scipy.constants import k, G, N_A
from matplotlib import pyplot as plt
plt.rc('text', usetex=True)

#We use Mars-like values from Yelle
#H in CO2 atmosphere, constant mixing ratio(p)

#initial composition control
K         = 3e6                                        #eddy diffusion coefficient, [cm^2/sec]
alpha     = -0.25                                      #from Yelle 
M_p       = 6.417e23                                   #mass of Mars, [kg]
R_p       = 3389.5e3                                   #radius of Mars, [m]
M_CO2     = 44.0095*1e-3                               #molar mass, [kg/mol]
M_H       = 1*1e-3                                     #molar mass, [kg/mol]
mole_frac = 1e-5                                       #mole fraction of H
m_a       = (mole_frac*M_H + (1-mole_frac)*M_CO2)/N_A  #mean molecular mass, [kg]
m_i       = M_H/N_A                                    #minor constituent mean molecular mass, [kg]

#initial grid control. the exobase temperature, T_inf, is left as an input parameter
z_0     = 80e3   #height above surface at bottom of the model, [m]
p_0     = 0.1    #pressure defining bottom of the model, [Pa]
p_inf   = 1e-6  #pressure defining roughly the exobase, [Pa]
p_steps = 10000   #atmosphere layers
T_0     = 100    #temperature at the bottom of the model (p_0), [K]

#retrieval of binary diffusion parameters
df = pd.read_csv("diffusion_coefficients.csv", comment="#")
#select a system here;
row = df[df["system"] == "H_in_CO2"].iloc[0]
A = row["A"]
s = row["s"]



#the implementation
def improved_limiting_flux(T_inf):
    '''calculates the improved limiting flux in [1/cm^2*s] for the above defined planet given the exobase temperature'''
    
    #1. define the p profile, derive T profile. use the Bates profile of Yelle eq. 19, but converted to function of pressure (not chi)
    p = np.logspace(np.log10(p_0), np.log10(p_inf), p_steps)           #pressure profile
    T = T_0 + (T_inf - T_0) * (1-((p/p_0)**0.75))                      #temperature profile
    
    #2. height profile using p and T(p) (integrating the HSE)
    z = np.zeros(p_steps)
    z[0] = z_0
    
    for i in range(0, p_steps-1):
        dp = p[i+1] - p[i]
        dz = -dp * ((k*T[i]*((z[i]+R_p)**2))/(m_a*p[i]*G*M_p))
        z[i+1] = z[i] + dz
    
    #3. number density profile, simply from ideal gass
    n = (p)/(k*T) * 1e-6  #[1/cm^3]
    
    #4. binary diffusion parameter and coefficient
    b = A*(T**s)  #[1/cm*s]
    D = b/n       #[cm^2/s]
    
    #5. g integral from Yelle
    #start with m_tilde needed for X_tilde
    
    T_grad  = - ((T[hom_arg+1]-T[hom_arg-1])/((p[hom_arg+1]-p[hom_arg-1]))) * p[hom_arg]
    m_tilde = m_i + (alpha*T_grad*(m_a/T[hom_arg]))