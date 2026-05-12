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
M_CO2     = 44.0095*1e-3                              #molar mass, [kg/mol]
M_H       = 1*1e-3                                    #molar mass, [kg/mol]
mole_frac = 1e-5                                      #mole fraction of H
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
def classical_limiting_flux(T_inf):
    '''calculates the classical limiting flux in [1/cm^2*s] for the above defined planet given the exobase temperature'''
    
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
    
    #4. homopause location
    b = A*(T**s)  #binary diffusion parameter, [1/cm*s]
    
    #now we locate b at the homopause by finding the altitude where D is closest to K
    D       = b/n                #[cm^2/s]
    diffs   = np.abs(K-D)
    hom_arg = np.argmin(diffs)
    b_hom   = b[hom_arg] * 1e2  #[1/m*s]
    
    #5. limiting flux calculation
    #start with scale height 
    H = (k*T[hom_arg]*((z[hom_arg]+R_p)**2))/(G*M_p*m_a)
    
    #then we need dT/d(xi) at the homopause, which via algebra is -dT/dp p
    T_grad  = - ((T[hom_arg+1]-T[hom_arg-1])/((p[hom_arg+1]-p[hom_arg-1]))) * p[hom_arg]
    m_tilde = m_i + (alpha*T_grad*(m_a/T[hom_arg]))
    
    flux_per_area = mole_frac * (b_hom/H) * (1-(m_tilde/m_a)) * 1e-4  #[1/cm^2*s]
    flux_tot      = flux_per_area * 4*np.pi*(((z[hom_arg]+R_p)*100)**2)      #[1/s]  !!!
    
    return(flux_per_area, flux_tot, n[-1])


#now we perform an exobase temperature sweep and plot the results
T_sweep = np.arange(100, 400)
flux_list = []
exobase_n_list = []

for temp in T_sweep:
    results = classical_limiting_flux(temp)
    flux_list.append(results[0])
    exobase_n_list.append(results[2])
    
fluxes = np.asarray(flux_list, dtype='float64')
exobase_n = np.asarray(exobase_n_list, dtype='float64')

plt.plot(T_sweep, fluxes)
plt.xlabel('Exobase temperature [K]')
plt.ylabel('Escape flux [cm$^{-2}$ $\cdot$ s$^{-1}$]')
plt.grid()
plt.savefig('plots/classical_lim.png')
plt.show()

#plt.plot(T_sweep, exobase_n)
#plt.xlabel('Exobase temperature [K]')
#plt.ylabel('Exobase number density [cm$^{-3}$]')
#plt.grid()
#plt.show()


    