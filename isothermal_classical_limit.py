import pandas as pd
import numpy as np
from scipy.constants import k, G, N_A
from scipy.integrate import quad


#Mars-like values from Jelle
#1% H in CO2 atmosphere, constant mixing ratio

T         = 100                                        #temperature profile, [K]
K         = 3e6                                        #eddy diffusion coefficient, [cm^2/sec]
M_p       = 6.417e23                                   #mass of Mars, [kg]
Mars_rad  = 3389.5e3                                   #radius of Mars, [m]
M_CO2     = 44.0095*10e-3                              #molar mass, [kg/mol]
M_H       = 1*10e-3                                    #molar mass, [kg/mol]
mole_frac = 10e-3                                      #mole fraction of H
m_a       = (mole_frac*M_H + (1-mole_frac)*M_CO2)/N_A  #mean molecular mass, [kg]
m_i       = M_H/N_A                                    #minor constituent mean molecular mass, [kg]



#1. binary diffusion parameter, b

df = pd.read_csv("diffusion_coefficients.csv", comment="#")
#select a system here;
row = df[df["system"] == "H_in_CO2"].iloc[0]
A = row["A"]
s = row["s"]
b = A*(T**s)  #[cgs]


#2. density profile and homopause location

#we follow Yelle and a Martian example and set the bottom of the atmosphere, z_0 at 80 km, with p_0 = 0.1 Pa
#the top of the atmosphere is set at 500 km

p_0   = 0.1                                               #pressure at bottom of the atmosphere [Pa]
z_0   = 80e3                                              #bottom of atmosphere [m]
z_top = 100e3                                             #top of atmosphere, [m]
steps = 1000                                              #layers of the atmosphere for computation 
z     = np.linspace(z_0, z_top, steps)                    #altitude grid, [m]


def density_profile(z):
    '''calculates the density profile of the atmosphere in particles/cm^3, with m_a, p_0, T, M_p and Mars_rad defined above''' #analytical solution
    rho_0 = (m_a*p_0)/(k*T)
    
    integrals = []
    for i in range(0, np.shape(z)[0]):
        I = ((G*M_p*m_a)/(k*T))*((1/(Mars_rad+z[i]))-(1/(Mars_rad+z[0]))) 
        integrals.append(I)
    
    ints = np.asarray(integrals, dtype='float64')
    return((rho_0 * np.exp(ints))/(m_a*10e6))


#now we locate the homopause by finding the altitude where D is closest to K
densities = density_profile(z)
D = b/densities
differences = np.abs(K-D)
homopause = z[np.argmin(differences)] #[m]
 
b_homopause = D[np.argmin(differences)]*10e2*densities[np.argmin(differences)]*10e6  #binary diffusion parameter at homopause, [1/m*s]


#3. limiting flux calculation

def classical_limiting_flux(b_hom, T_hom, z_hom):
    '''calculates the classical limiting flux of the species at question for an isothermal atmosphere'''
    H_a = (k*T_hom*(z_hom**2))/(G*M_p*m_a)   #the scale height at the homopause, [m]
    
    return ((mole_frac*(b_hom/H_a)*(1-(m_i/m_a))), (mole_frac*(b_hom/H_a)*(1-(m_i/m_a)) * 4*np.pi*((z_hom+Mars_rad)**2)))  #flux in particles/second*m^2 and flux from entire planet in particles/second

res = classical_limiting_flux(b_homopause, T, homopause)
print(f"the classical limiting flux is: {res[0]:g} particles per second per square meter, or {res[1]:g} particles per second overall")


    
