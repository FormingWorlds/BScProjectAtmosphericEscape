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
    
    #5. we need to find the homopause and the exobase to be used as our integral limits
    #exobase
    H = (k*T*((R_p+z)**2))/(m_a*G*M_p)   #scale height, [m]
    mfp = (D*1e-4) * (m_a/(k*T))**0.5    #mean free path, [m] 
    exo_diffs = np.abs(H-mfp)
    exo_arg = np.argmin(exo_diffs)
   
    #homopause
    hom_diffs   = np.abs(K-D)
    hom_arg = np.argmin(hom_diffs)
    z_hom = z[hom_arg]
    
    #lets re-define the arrays with the new limits for easier use
    T = T[hom_arg:exo_arg+1]
    p = p[hom_arg:exo_arg+1]
    xi = -np.log(p/p[0])        #the definition of xi from Yelle
    
    #and the coefficients in SI
    D = D*1e-4
    K_si = K*1e-4
    
    #6. g integral from Yelle
    #can do all in a single loop: need m_tilde for X_tilde integral, which is needed for g integral
    int_1 = 0   #exponential integral in X_i expression
    g_int = 0   #g function integral
    
    for i in range(0, np.shape(T)[0]-1):
        #differential in xi(i)
        d_xi = xi[i+1]-xi[i]
        
        #m_tilde(i) calculation
        T_grad  = - ((T[i+1]-T[i])/((p[i+1]-p[i]))) * p[i]
        m_tilde = m_i + (alpha*T_grad*(m_a/T[i]))
        
        #X_i_tilde(i) calculation
        int_1 = int_1 + (((1-(m_tilde/m_a)) * (D[i]/(D[i]+K_si))) * d_xi)
        X_i_tilde = mole_frac * np.exp(int_1)
        
        #g(i) calculation
        r_0 = (R_p+z_hom)**2    #homopause radius squared
        g_int = g_int + ((k*T[i]*r_0)/(X_i_tilde*(n[i]*1e6)*G*M_p*m_a*(D[i]+K_si))) * d_xi
        
    #7. final escape flux, inverse of g
    flux = (1/g_int)*1e-4     #[1/cm^s*sec]
    
    return(flux, T[-1])


#now we perform an exobase temperature sweep and plot the results
T_sweep = np.arange(100, 400)
flux_list = []
T_list = []

for temp in T_sweep:
    results = improved_limiting_flux(temp)
    flux_list.append(results[0])
    T_list.append(results[1])
    
fluxes = np.asarray(flux_list, dtype='float64')
temps = np.asarray(T_list, dtype='float64')

plt.plot(temps, fluxes)
plt.xlabel('Exobase temperature [K]')
plt.ylabel('Escape flux [cm$^{-2}$ $\cdot$ s$^{-1}$]')
plt.grid()
plt.savefig('plots/improved_lim.png')
plt.show()