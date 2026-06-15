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
            
            
            
            p_b_list.append(p_b)
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
            
        T_b_arr = np.asarray(T_b_list, dtype='float64')
        p_b_arr = np.asarray(p_b_list, dtype='float64')
        z_b_arr = np.asarray(z_b_list, dtype='float64')
        rho_b_arr = np.asarray(rho_b_list, dtype='float64')
        exo_id_arr = np.asarray(exo_id_list, dtype='float64')
        H_arr = np.asarray(H_list, dtype='float64')
        mfp_arr = np.asarray(mfp_list, dtype='float64')
        return(p_b_arr, T_b_arr, T_inf, hom_list, exo_id_arr, H_arr, mfp_arr)
            


        #     if isinstance(exo_id, np.int64):
        #         #limiting flux calculation
        #         #minor constituent (here - atomic hydrogen) mean molecular mass, [kg]
        #         m_i = molarmasses['H']/N_A
                
        #         #re-defining the arrays with the new limits for easier use
        #         z_hom = z_b[hom_id]
        #         z_exo = z_b[exo_id+1]
        #         T_b = T_b[hom_id:exo_id+1]
        #         p_b = p_b[hom_id:exo_id+1]
        #         D_b = D[hom_id:exo_id+1]
        #         rho_b = rho_b[hom_id:exo_id+1]
        #         K_b = K_b[hom_id:exo_id+1]
        #         mmw_b = mmw_b[hom_id:exo_id+1]
        #         for specie in VMR_b:
        #             VMR_b[specie] = VMR_b[specie][hom_id:exo_id+1]
                    
        #         #extending D
        #         p_atm = p_b/101325 #pressure in [atm], needed for Slattery function
        #         m_a = mmw_b/N_A #overall mean molecular mass of a single molecule, [kg]
        #         n_b = (rho_b/m_a) * 1e-4  #particle number density in cgs [1/cm^3]
                
        #         A = atm.diff_coeffs['H']['A']
        #         s = atm.diff_coeffs['H']['s']
                
        #         if np.shape(atm.diff_coeffs['H']['A'])[0] == 0:
        #             D_b = slattery_diff(T_b, p_atm, 'H', system)
        #         else:
        #             b = A*(T_b**s)
        #             D_b = b/n_b
                    
                    
        #         D_list.append(D_b)
        #         p_b_list.append(p_b)
                    
                
        #         #the definition of xi from Yelle
        #         xi = -np.log(p_b/p_b[0])
                
        #         #re-defining the coefficients in SI
        #         D_b = D_b*1e-4
        #         K_si = K_b*1e-4
                
                
        #         #g integral from Yelle, now over all species taking into account their dissociation fraction
        #         flux_from_specie = []
        #         for H_specie in disso_fracs:
        #             if(VMR_b[H_specie][0] > 1e-5):
        #                 int_1 = 0   #exponential integral in X_i expression
        #                 g_int = 0   #g function integral
        #                 mole_frac_homo = ((VMR_b[H_specie][0] * disso_fracs[H_specie] * Natoms[H_specie] * Hfrac[H_specie])/((1 - (VMR_b[H_specie][0] * disso_fracs[H_specie]))+(VMR_b[H_specie][0] * disso_fracs[H_specie] * Natoms[H_specie])))
        #                 for i in range(0, np.shape(T_b)[0]-1):
        #                     #differential in xi(i)
        #                     d_xi = xi[i+1]-xi[i]
                            
        #                     #m_tilde(i) calculation
        #                     T_grad  = (T_b[i+1]-T_b[i])/d_xi
        #                     m_tilde = m_i + (alpha*T_grad*(m_a/T_b[i]))
                            
        #                     #X_i_tilde(i) calculation
        #                     int_1 = int_1 + (((1-(m_tilde/m_a)) * (D_b[i]/(D_b[i]+K_si))) * d_xi)
        #                     X_i_tilde = mole_frac_homo * np.exp(int_1)
                            
        #                     #g(i) calculation
        #                     r_0 = (R+z_hom)**2    #homopause radius squared
        #                     g_int = g_int + ((k*T_b[i]*r_0)/(X_i_tilde*(n_b[i]*1e4)*G*M*m_a*(D_b[i]+K_si))) * d_xi

                        
        #                 #final escape flux, inverse of g
        #                 flux_from_specie.append((1/g_int)*1e-4)
                    

            
                
                
                
        #         flux_from_specie_arr = np.asarray(flux_from_specie, dtype='float64')
        #         flux = np.sum(flux_from_specie_arr)
        #         flux_list.append(flux)
        #         #flux_arr = np.asarray(flux_list, dtype='float64')
                
        #         flux_SI_list.append(flux * 10000 * (4*np.pi*((R+z_exo)**2)) * m_i)
        #         #flux_SI_arr = np.asarray(flux_SI_list, dtype='float64')
                
        #         exo_height_list.append(z_exo/1000)
        #         exo_pressure_list.append(p_b[-1]/100000)
        #         exo_temp_list.append(T_b[-1])
                
        #         hom_height_list.append(z_hom/1000)
        #         hom_pressure_list.append(p_b[0]/100000)
        #         hom_temp_list.append(T_b[0])

        
        # flux_arr = np.asarray(flux_list, dtype='float64')
        # flux_SI_arr = np.reshape((np.asarray(flux_SI_list, dtype='float64')), (1, np.shape(T_inf)[0]))
              
        # exo_height_arr = np.asarray(exo_height_list, dtype='float64')
        # exo_pressure_arr = np.asarray(exo_pressure_list, dtype='float64')
        # exo_temp_arr = np.asarray(exo_temp_list, dtype='float64')
        
        # hom_height_arr = np.asarray(hom_height_list, dtype='float64')
        # hom_pressure_arr = np.asarray(hom_pressure_list, dtype='float64')
        # hom_temp_arr = np.asarray(hom_temp_list, dtype='float64')
        
        
        
            
        # #return(D_list, p_b_list, T_inf)    
        # return(T_inf, flux_arr, flux_SI_arr, hom_height_arr, exo_height_arr, hom_pressure_arr, exo_pressure_arr, hom_temp_arr, exo_temp_arr)