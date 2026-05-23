import numpy as np
from scipy.constants import k, G, N_A
from slattery_diff import slattery_diff


alpha     = -0.25       #from Yelle 


#the function needs to be able to do the following:
#1. each of the 6 diffusing species has a changeable "dissociation" parameter controlling the amount of free hydrogen diffusing
#2. the diffusion coefficicent of each species can be changed (we will initially assume H-systtem for all species, since theyre dissociated)


def proteus_improved_limiting_flux(atm, system, disso_fracs):
        """atm = class
        system = string of main consituent
        disso_fracs = library of disso fracs"""
        
        
        T = atm.temperature
        
        p = atm.pressure
        K = atm.Kzz
        
        rho = atm.density
        mmw = atm.mmw
        z = atm.height
        
        
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
            
        
            
        H = atm.scale_height()        
        mfp = (D*1e-4) * ((m_a/(k*T))**0.5)    #mean free path, [m] 
        
        
        Kres = atm.Kzz_extension()
        
        Tres = atm.isothermal_extension()
            
            
        return(Tres[0], Tres[1], Tres[2], Tres[3], Tres[4], Kres[0], Kres[1])