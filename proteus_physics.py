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
        
        A = atm.diff_coeffs['H']['A']
        s = atm.diff_coeffs['H']['s']
        T = atm.temperature
        p = atm.pressure
        p_atm = p/101325 #[atm]
        n = (atm.density / (atm.mmw/N_A)) * 1e-4
        
        b = A * (T**(s))
        D = b/n
        
        
        slattery = slattery_diff(T, p_atm, 'H', system)
        
        
        return(D, slattery)
