import numpy as np
from scipy.constants import k, G, N_A
from slattery_diff import slattery_diff
from proteus_extra_functions import H_mean_free_path, H_scale_height


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
        
        return(mfp, H, p, D, K)
        
        
        # #to go on im just gonna pretend as if i have the homopause and exobase indices, hom_id and exo_id
        # exo_arg = np.argmin(exo_diffs)
        # hom_arg = np.argmin(hom_diffs)
        
        # #minor constituent (here - atomic hydrogen) mean molecular mass, [kg]
        # m_i = molarmasses['H']/N_A
        
        # #re-defining the arrays with the new limits for easier use
        # z_hom = z[hom_id]
        # T = T[hom_id:exo_id+1]
        # p = p[hom_id:exo_id+1]
        # D = D[hom_id:exo_id+1]
        # n = n[hom_id:exo_id+1]
        
        # for specie in VMR:
        #     VMR[specie] = VMR[specie][hom_id:exo_id+1]
        
        # #the definition of xi from Yelle
        # xi = -np.log(p/p[0])
        
        # #re-defining the coefficients in SI
        # D = D*1e-4
        # K_si = K*1e-4
        
        # #g integral from Yelle, now over all species taking into account their dissociation fraction
        # g_list = []
        # for H_specie in disso_fracs:
        #     int_1 = 0   #exponential integral in X_i expression
        #     g_int = 0   #g function integral
        #     mole_frac_homo = ((VMR[H_specie][0] * disso_fracs[H_specie] * Natoms[H_specie] * Hfrac[H_specie])/((1 - (VMR[H_specie][0] * disso_fracs[H_specie]))+(VMR[H_specie][0] * disso_fracs[H_specie] * Natoms[H_specie])))
        #     for i in range(0, np.shape(T)[0]-1):
        #         #differential in xi(i)
        #         d_xi = xi[i+1]-xi[i]
                
        #         #m_tilde(i) calculation
        #         T_grad  = - ((T[i+1]-T[i])/((p[i+1]-p[i]))) * p[i]
        #         m_tilde = m_i + (alpha*T_grad*(m_a/T[i]))
                
        #         #X_i_tilde(i) calculation
        #         int_1 = int_1 + (((1-(m_tilde/m_a)) * (D[i]/(D[i]+K_si))) * d_xi)
        #         X_i_tilde = mole_frac_homo * np.exp(int_1)
                
        #         #g(i) calculation
        #         r_0 = (R+z_hom)**2    #homopause radius squared
        #         g_int = g_int + ((k*T[i]*r_0)/(X_i_tilde*(n[i]*1e4)*G*M*m_a*(D[i]+K_si))) * d_xi
        #     g_list.append(g_int)
        
        # g_array = np.asarray(g_list, dtype='float64')
        # g_tot = np.sum(g_array)
        
        
        # #final escape flux, inverse of g
        # flux = (1/g_tot)*1e-4     #[1/cm^s*sec]
        
        # return(flux)