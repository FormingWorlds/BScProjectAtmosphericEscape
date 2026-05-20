import pandas as pd
import numpy as np
from scipy.constants import k, G, N_A
from slattery_diff import slattery_diff



def improved_limiting_flux(T_inf, minor_const, major_const, M_minor, M_major, mole_frac_homop, K, M_p, R_p, z_0, p_0, p_inf, p_steps, T_0, alpha):
    '''calculates the improved limiting flux in [1/cm^2*s] and exobase temperature [K] for the above defined planet given the top layer temperature (approx exobase), major and 
    minor constituents (str) and their molar masses [kg/mol], the minor constituent mole fraction at homopause and a constant eddy diffusion parameter [cm^2/s]'''
    
        
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
    
    #5. binary diffusion parameter and coefficient
    df = pd.read_csv("diffusion_coefficients.csv")
    row = df[(df["minor_const"] == minor_const) & (df["major_const"] == major_const)]
    A = row["A"].values
    s = row["s"].values
        
    if np.shape(A)[0] == 0:
        p_atm = p/101325 #[atm]
        D = slattery_diff(T, p_atm, minor_const, major_const)
    else:
        b = A*(T**s)  #[1/cm*s]
        D = b/n       #[cm^2/s]
    
    #6. we need to find the homopause and the exobase to be used as our integral limits
    #exobase
    H = (k*T*((R_p+z)**2))/(m_a*G*M_p)   #scale height, [m]
    mfp = (D*1e-4) * (m_a/(k*T))**0.5    #mean free path, [m] 
    exo_diffs = np.abs(H-mfp)
    exo_arg = np.argmin(exo_diffs)
   
    #homopause
    hom_diffs   = np.abs(K-D)
    hom_arg = np.argmin(hom_diffs)
    z_hom = z[hom_arg]
    
    #check if exobase and homopause were found
    if(exo_arg <= hom_arg):
        return f"ERROR: homopause index: {hom_arg}, Kmax: {K:g}, exobase index: {exo_arg}, Dmax: {D[-1]:g}, H_0: {H[0]:g}, H_max: {H[-1]:g}, mfp_0: {mfp[0]:g}, mfp_max: {mfp[-1]:g}"
    
    #lets re-define the arrays with the new limits for easier use
    T = T[hom_arg:exo_arg+1]
    p = p[hom_arg:exo_arg+1]
    D = D[hom_arg:exo_arg+1]
    n = n[hom_arg:exo_arg+1]
    z = z[hom_arg:exo_arg+1]
    
    xi = -np.log(p/p[0])        #the definition of xi from Yelle
    
    #and the coefficients in SI
    D = D*1e-4
    K_si = K*1e-4
    
    #7. g integral from Yelle
    #can do all in a single loop: need m_tilde for X_tilde integral, which is needed for g integral
    int_1 = 0   #exponential integral in X_i expression
    g_int = 0   #g function integral
    
    # X_i_list = []
    # g_list = []
    
    for i in range(0, np.shape(T)[0]-1):
        #differential in xi(i)
        d_xi = xi[i+1]-xi[i]
        
        #m_tilde(i) calculation
        T_grad  = - ((T[i+1]-T[i])/((p[i+1]-p[i]))) * p[i]
        m_tilde = m_i + (alpha*T_grad*(m_a/T[i]))
        
        #X_i_tilde(i) calculation
        int_1 = int_1 + (((1-(m_tilde/m_a)) * (D[i]/(D[i]+K_si))) * d_xi)
        X_i_tilde = mole_frac_homop * np.exp(int_1)
        
        
        #g(i) calculation
        r_0 = (R_p+z_hom)**2    #homopause radius squared
        g_int = g_int + ((k*T[i]*r_0)/(X_i_tilde*(n[i]*1e6)*G*M_p*m_a*(D[i]+K_si))) * d_xi
        
        
        # X_i_list.append(X_i_tilde)
        # g_list.append(g_int)
        
    # X_i_exo = np.asarray(X_i_list, dtype='float64')
    # g_exo = np.asarray(g_list, dtype='float64')
    
    # X_i_func = X_i_exo * (1-(g_exo*g_int))
    
        
    #8. final escape flux, inverse of g
    flux = (1/g_int)*1e-4     #[1/cm^s*sec]
    
    return(flux, T[-1], hom_arg, exo_arg, D[-1], H[0], H[-1], mfp[0], mfp[-1])