import numpy as np
from scipy.constants import k, G, N_A
from slattery_diff import slattery_diff
from proteus_extra_functions import scale_height


alpha     = -0.25       #from Yelle 


#the function needs to be able to do the following:
#1. each of the 6 diffusing species has a changeable "dissociation" parameter controlling the amount of free hydrogen diffusing
#2. the diffusion coefficicent of each species can be changed (we will initially assume H-systtem for all species, since theyre dissociated)


def proteus_improved_limiting_flux(atm, system, disso_fracs):
        """atm = class
        system = string of main consituent
        disso_fracs = library of disso fracs"""
        
        Kres = atm.Kzz_extension()
        Tres = atm.isothermal_extension()
        
        T = Tres[0]
        VMR = Tres[1]
        mmw = Tres[2]
        rho = Tres[3]
        z = Tres[4]
        
        K = Kres[0]
        p = Kres[1]
        
        M = atm.planet_mass
        R = atm.planet_rad
        
        
        m_a = mmw/N_A
        p_atm = p/101325 #[atm]
        n = (rho/m_a) * 1e-4
        
        
        A = atm.diff_coeffs['H']['A']
        s = atm.diff_coeffs['H']['s']
        
        if np.shape(atm.diff_coeffs['H']['A'])[0] == 0:
            D = slattery_diff(T, p_atm, 'H', system)
        else:
            b = A*(T**s)
            D = b/n
            
        
            
             
        mfp = (D*1e-4) * ((m_a/(k*T))**0.5)    #mean free path, [m] 
        
        H = scale_height(T, R, z, M, mmw) 
            
            
        return(T, VMR, mmw, rho, z, K, p, H, mfp, D)