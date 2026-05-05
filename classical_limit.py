import pandas as pd
import numpy as np
from scipy.constants import k, G, N_A

#assumes constant mole fraction

#1. composition control

#We initially use Mars-like values from Yelle
#H in CO2 atmosphere, constant mixing ratio(p)
K         = 3e6                                        #eddy diffusion coefficient, [cm^2/sec]
alpha     = -0.25                                      #from Yelle 
M_p       = 6.417e23                                   #mass of Mars, [kg]
R_p       = 3389.5e3                                   #radius of Mars, [m]
M_CO2     = 44.0095*10e-3                              #molar mass, [kg/mol]
M_H       = 1*10e-3                                    #molar mass, [kg/mol]
mole_frac = 10e-5                                      #mole fraction of H
m_a       = (mole_frac*M_H + (1-mole_frac)*M_CO2)/N_A  #mean molecular mass, [kg]
m_i       = M_H/N_A                                    #minor constituent mean molecular mass, [kg]


#2. grid control.  
#we define the p profile, derive T profile. use the Bates profile of Yelle eq. 19, but converted to function of pressure (not chi)

z_0     = 80e3   #height above surface at bottom of the model, [m]
p_0     = 0.1    #pressure defining bottom of the model, [Pa]
p_inf   = 10e-6  #pressure defining roughly the exobase, [Pa]
p_steps = 10000   #atmosphere layers
T_0     = 100    #temperature at the bottom of the model (p_0), [K]
T_inf   = 250    #exobase temperature (p_inf), [K]

p = np.linspace(p_0, p_inf, p_steps)           #pressure profile
T = T_0 + (T_inf - T_0) * (1-((p/p_0)**0.75))  #temperature profile


#3. height and number density profiles using p and T(p) (and potentially m(p))

def height_profile(p, T, m, M_p):
    '''gives the corresponding height profile given the p, T and m profiles and the planet mass, in meters'''
    z = np.zeros(p_steps)
    z[0] = z_0
    
    for i in range(0, p_steps-1):
        dp = p[i+1] - p[i]
        dz = -dp * ((k*T[i]*((z[i]+R_p)**2))/(m*p[i]*G*M_p))
        z[i+1] = z[i] + dz
        
    return z

z = height_profile(p, T, m_a, M_p)

#number density profile, simply from ideal gass
n = (p)/(k*T) * 10e-6  #[1/cm^3]


#4. homopause location
#retrieval of binary diffusion parameters

df = pd.read_csv("diffusion_coefficients.csv", comment="#")
#select a system here;
row = df[df["system"] == "H_in_CO2"].iloc[0]
A = row["A"]
s = row["s"]
b = A*(T**s)  #[1/cm*s]

#now we locate b at the homopause by finding the altitude where D is closest to K
D = b/n       #[cm^2/s]
diffs = np.abs(K-D)
hom_arg = np.argmin(diffs)
b_hom = b[hom_arg] * 10e2  #[1/m*s]


#5. limiting flux calculation

def classical_limiting_flux(X_i, b_hom, T, z, M_p, m_a, m_i):
    '''finds the classical limiting flux'''
    #start with scale height 
    H = (k*T[0]*((z[0]+R_p)**2))/(G*M_p*m_a)
    
    #then we need dT/d(xi) at the homopause, which turns out is -dT/dp p
    T_grad = - ((T[hom_arg+1]-T[hom_arg-1])/((p[hom_arg+1]-p[hom_arg-1]))) * p[hom_arg]
    m_tilde = m_i + (alpha*T_grad*(m_a/T[hom_arg]))
    
    flux_per_area = X_i * (b_hom/H) * (1-(m_tilde/m_a)) * 10e-4   #[1/cm^2*s]
    flux_tot = flux_per_area * 4*np.pi*((z[hom_arg]+R_p)**2)      #[1/s]
    return(flux_per_area, flux_tot)

results = classical_limiting_flux(mole_frac, b_hom, T, z, M_p, m_a, m_i)
print(f'flux = {results[0]:g} [1/cm^2 * s]')


    