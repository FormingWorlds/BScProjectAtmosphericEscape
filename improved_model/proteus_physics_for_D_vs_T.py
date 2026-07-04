import numpy as np
from scipy.constants import k, G, N_A
from slattery_diff import slattery_diff
from proteus_extra_functions import H_mean_free_path, H_scale_height, homopause_index_finder, Bates_extension, exobase_index_finder


alpha     = -0.25       #from Yelle 


#each of the 6 diffusing species has a changeable "dissociation" parameter controlling the amount of free hydrogen diffusing
#fraction of the initial molecules that are fully dissociated (the physics module takes into account the differing amount of atoms in different molecules)
disso_fracs = {   
    "H2" : 1,
    "H2O": 1,
    "H"  : 1,
    "CH4": 1,
    "NH3": 1,
    "H2S": 1
}


def proteus_profiles(atm, system, disso_fracs):
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
    hom_list = []
        
    #if homopause found - continue, if not returns None
    if isinstance(hom_id, np.int64):
        
        flux_list = []
        flux_SI_list = []
        exo_height_list = []
        exo_pressure_list = []
        exo_temp_list = []
        hom_height_list = []
        hom_pressure_list = []
        hom_temp_list = []
        
        T_inf = np.array([200, 300, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000]) #[K]
        
        #atmosphere profile plotting outputs
        p_b_list = []
        T_b_list = []
        z_b_list = []
        rho_b_list = []
        exo_id_list = []
        
        mfp_list = []
        H_list = []
        
        D_list = []
        
        hom_list.append(hom_id)
        
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
            
            
            
            #p_b_list.append(p_b)
            T_b_list.append(T_b)
            z_b_list.append(z_b)
            rho_b_list.append(rho_b)
            

            
            #now we find the exobase index, where mfp >= H
            mfp = H_mean_free_path(VMR_b, rho_b, mmw_b)
            H = H_scale_height(T_b, R, z_b, M, molarmasses)
            exo_id = exobase_index_finder(H, mfp)
            
            mfp_list.append(mfp)
            H_list.append(H)
            exo_id_list.append(exo_id)

            if isinstance(exo_id, np.int64):
                #limiting flux calculation
                #minor constituent (here - atomic hydrogen) mean molecular mass, [kg]
                m_i = molarmasses['H']/N_A
                
                #re-defining the arrays with the new limits for easier use
                z_hom = z_b[hom_id]
                z_exo = z_b[exo_id+1]
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
                    
                    
                D_list.append(D_b)
                p_b_list.append(p_b)
                
        return(D_list, p_b_list)