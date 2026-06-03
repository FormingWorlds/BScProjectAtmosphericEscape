import numpy as np
from scipy.constants import k, G, N_A
from slattery_diff import slattery_diff
from proteus_extra_functions import H_mean_free_path, H_scale_height, homopause_index_finder, Bates_extension, exobase_index_finder


alpha     = -0.25       #from Yelle 


#the function needs to be able to do the following:
#each of the 6 diffusing species has a changeable "dissociation" parameter controlling the amount of free hydrogen diffusing


def proteus_improved_limiting_flux(atm, system, disso_fracs):
    """atm = class
    system = string of main consituent
    disso_fracs = library of disso fracs"""
    
    #calling for all needed values
    M = atm.planet_mass
    R = atm.planet_rad
    molarmasses = atm.molar_masses
    Natoms = atm.Natoms
    Hfrac = atm.Hfrac
    
    Kres = atm.Kzz_extension()
    Tres = atm.isothermal_extension()
    
    T = Tres[0]
    VMR = Tres[1]
    mmw = Tres[2]
    rho = Tres[3]
    z = Tres[4]
    
    K = Kres[0]
    p = Kres[1]
    
    #some definitions
    m_a = mmw/N_A #overall mean molecular mass of a single molecule, [kg]
    p_atm = p/101325 #pressure in [atm], needed for Slattery function
    n = (rho/m_a) * 1e-4  #particle number density in cgs [1/cm^3]
    
    #diffusion parameter fetch for homopause location and limit calculation
    A = atm.diff_coeffs['H']['A']
    s = atm.diff_coeffs['H']['s']
    
    if np.shape(atm.diff_coeffs['H']['A'])[0] == 0:
        D = slattery_diff(T, p_atm, 'H', system)
    else:
        b = A*(T**s)
        D = b/n
        
    #hydrogen mfp and scale height for exobase location  
    mfp = H_mean_free_path(VMR, rho, mmw)
    H = H_scale_height(T, R, z, M, molarmasses)
    
    #homopause location
    hom_id = homopause_index_finder(K, D, H, mfp)
        
    #if homopause found - continue, if not returns None
    if isinstance(hom_id, np.int64):
        flux_list = []
        T_inf = np.array([100, 300, 500, 1000, 2000, 3000]) #[K]
        for temp in T_inf:               
            #extending as a Bates profile above homopause 
            Bates = Bates_extension(hom_id, temp, M, R, T, p, z, rho, mmw, VMR, K)
            T_b = Bates[0]
            p_b = Bates[1]
            z_b = Bates[2]
            rho_b = Bates[3]
            mmw_b = Bates[4]
            VMR_b = Bates[5]
            K_b = Bates[6]
            
            #now we find the exobase index, where mfp >= H
            mfp = H_mean_free_path(VMR_b, rho_b, mmw_b)
            H = H_scale_height(T_b, R, z_b, M, molarmasses)
            exo_id = exobase_index_finder(H, mfp)
            
            
            #limiting flux calculation
            #minor constituent (here - atomic hydrogen) mean molecular mass, [kg]
            m_i = molarmasses['H']/N_A
            
            #re-defining the arrays with the new limits for easier use
            z_hom = z_b[hom_id]
            T_b = T_b[hom_id:exo_id+1]
            p_b = p_b[hom_id:exo_id+1]
            D_b = D[hom_id:exo_id+1]
            rho_b = rho_b[hom_id:exo_id+1]
            K_b = K_b[hom_id:exo_id+1]
            mmw_b = mmw_b[hom_id:exo_id+1]
            for specie in VMR_b:
                VMR_b[specie] = VMR_b[specie][hom_id:exo_id+1]
                
            #extending D
            p_atm = p_b/101325 #pressure in [atm], needed for Slattery function
            m_a = mmw_b/N_A #overall mean molecular mass of a single molecule, [kg]
            n_b = (rho_b/m_a) * 1e-4  #particle number density in cgs [1/cm^3]
            
            A = atm.diff_coeffs['H']['A']
            s = atm.diff_coeffs['H']['s']
            
            if np.shape(atm.diff_coeffs['H']['A'])[0] == 0:
                D_b = slattery_diff(T_b, p_atm, 'H', system)
            else:
                b = A*(T_b**s)
                D_b = b/n_b
            
            #the definition of xi from Yelle
            xi = -np.log(p_b/p_b[0])
            
            #re-defining the coefficients in SI
            D_b = D_b*1e-4
            K_si = K_b*1e-4
            
            #g integral from Yelle, now over all species taking into account their dissociation fraction
            flux_from_specie = []
            for H_specie in disso_fracs:
                if(disso_fracs[H_specie] > 1e-5):
                    int_1 = 0   #exponential integral in X_i expression
                    g_int = 0   #g function integral
                    mole_frac_homo = ((VMR_b[H_specie][0] * disso_fracs[H_specie] * Natoms[H_specie] * Hfrac[H_specie])/((1 - (VMR_b[H_specie][0] * disso_fracs[H_specie]))+(VMR_b[H_specie][0] * disso_fracs[H_specie] * Natoms[H_specie])))
                    for i in range(0, np.shape(T_b)[0]-1):
                        #differential in xi(i)
                        d_xi = xi[i+1]-xi[i]
                        
                        #m_tilde(i) calculation
                        T_grad  = (T_b[i+1]-T_b[i])/d_xi
                        m_tilde = m_i + (alpha*T_grad*(m_a/T_b[i]))
                        
                        #X_i_tilde(i) calculation
                        int_1 = int_1 + (((1-(m_tilde/m_a)) * (D_b[i]/(D_b[i]+K_si))) * d_xi)
                        X_i_tilde = mole_frac_homo * np.exp(int_1)
                        
                        #g(i) calculation
                        r_0 = (R+z_hom)**2    #homopause radius squared
                        g_int = g_int + ((k*T_b[i]*r_0)/(X_i_tilde*(n_b[i]*1e4)*G*M*m_a*(D_b[i]+K_si))) * d_xi
                    
                    #final escape flux, inverse of g
                    flux_from_specie.append((1/g_int)*1e-4)
                
            flux_from_specie_arr = np.asarray(flux_from_specie, dtype='float64')
            flux = np.sum(flux_from_specie_arr)
            flux_list.append(flux)
            flux_arr = np.asarray(flux_list, dtype='float64')
                
        return(flux_arr, T_inf)