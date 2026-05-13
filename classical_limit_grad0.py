import pandas as pd
import numpy as np
from scipy.constants import k, G, N_A


#planetary bulk parameter control
alpha     = -0.25                                      #from Yelle 
M_p       = 6.417e23                                   #mass of Mars, [kg]
R_p       = 3389.5e3                                   #radius of Mars, [m]

#initial grid control. the temperature at the top layer, T_inf, is left as an input parameter
z_0     = 80e3   #height above surface at bottom of the model, [m]
p_0     = 0.1    #pressure defining bottom of the model, [Pa]
p_inf   = 1e-6  #pressure defining roughly the exobase, [Pa]
p_steps = 10000   #atmosphere layers
T_0     = 100    #temperature at the bottom of the model (p_0), [K]



def classical_limiting_flux_grad0(T_inf, minor_const, major_const, M_minor, M_major, mole_frac_homop, K):
    '''calculates the classical limiting flux in [1/cm^2*s] for the above defined planet given the top layer temperature (approx exobase), major and 
    minor constituents (str) and their molar masses [kg/mol], the minor constituent mole fraction at homopause and a constant eddy diffusion parameter [cm^2/s], 
    assuming the temperature gradient is 0 (negligible)'''
    
    #0. retrieval of binary diffusion parameters
    df = pd.read_csv("diffusion_coefficients.csv")
    row = df[(df["minor_const"] == minor_const) & (df["major_const"] == major_const)]
    A = row["A"].values
    s = row["s"].values
    
    #1. molecular masses
    m_a = (mole_frac_homop*M_minor + (1-mole_frac_homop)*M_major)/N_A  #overall mean molecular mass, [kg]
    m_i = M_minor/N_A                                      #minor constituent mean molecular mass, [kg]
    
    #2. define the p profile, derive T profile. use the Bates profile of Yelle eq. 19, but converted to function of pressure (not chi)
    p = np.logspace(np.log10(p_0), np.log10(p_inf), p_steps)           #pressure profile
    T = T_0 + (T_inf - T_0) * (1-((p/p_0)**0.75))                      #temperature profile
    
    #3. height profile using p and T(p) (integrating the HSE)
    z = np.zeros(p_steps)
    z[0] = z_0
    
    for i in range(0, p_steps-1):
        dp = p[i+1] - p[i]
        dz = -dp * ((k*T[i]*((z[i]+R_p)**2))/(m_a*p[i]*G*M_p))
        z[i+1] = z[i] + dz
    
    #4. number density profile, simply from ideal gass
    n = (p)/(k*T) * 1e-6  #[1/cm^3]
    
    #5. homopause location
    b = A*(T**s)  #binary diffusion parameter, [1/cm*s]
    
    #now we locate b at the homopause by finding the altitude where D is closest to K
    D       = b/n                #[cm^2/s]
    diffs   = np.abs(K-D)
    hom_arg = np.argmin(diffs)
    b_hom   = b[hom_arg] * 1e2  #[1/m*s]
    
    #5. limiting flux calculation
    #start with scale height 
    H = (k*T[hom_arg]*((z[hom_arg]+R_p)**2))/(G*M_p*m_a)
    
    m_tilde = m_i
    
    flux = mole_frac_homop * (b_hom/H) * (1-(m_tilde/m_a)) * 1e-4  #[1/cm^2*s]
    
    return(flux)    