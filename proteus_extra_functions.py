import numpy as np
from scipy.constants import k, G, N_A


def H_scale_height(T, R, z, M, molar_masses):
    """finds atomic hydrogen density scale height given T grid, height grid, planet mass and radius and H molar mass"""
    
    H = (k * T * ((R + z)**2))/(G * M * (molar_masses["H"]/N_A))
    return(H)


kinetic_radii = {
    "H2" : 289/2*1e-12, #kinetic D
    "H2O": 265/2*1e-12, #kinetic D
    "CH4": 380/2*1e-12, #kinetic D
    "NH3": 260/2*1e-12, #kinetic D
    "H2S": 360/2*1e-12, #kinetic D
    "N2": 364/2*1e-12, #kinetic D
    "CO2": 330/2*1e-12, #kinetic D
    "CO" : 376/2*1e-12, #kinetic D
    "SO2" : 360/2*1e-12, #kinetic D
    "O2" : 346/2*1e-12, #kinetic D
    "O" : 304/2*1e-12, #vanderwaals d
    "S" : 360/2*1e-12, #vanderwaals d
    "N" : 310/2*1e-12, #vanderwaals d
    "H"  : 218/2*1e-12 #vanderwaals d
}


def H_mean_free_path(vmr_dict, rho, mmw):
    """finds the atomic hydrogen mean free path in meters given mole fraction of species, density and mmw"""
    
    mfp_list = []
    for i in range(0, np.size(rho)):
        this_inv_mfp = 0
        for specie in vmr_dict:
            if(vmr_dict[specie][i] > 0.005):
                sigma = np.pi * ((kinetic_radii["H"] + kinetic_radii[specie])**2)
                specie_density = (rho[i] * vmr_dict[specie][i])/(mmw[i]/N_A)
                this_inv_mfp = this_inv_mfp + (sigma*specie_density)
        mfp_list.append(1/this_inv_mfp)
    
    mfp = np.asarray(mfp_list, dtype='float64')
    return(mfp)


def homopause_index_finder(K, D, H, mfp):
    """finds and returns the homopause index (where D becomes larger than K) if it exists"""

    if(np.shape(np.where(D >= K)[0])[0] == 0):
        hom_id = 'homopause not found'
    elif((np.shape(np.where(mfp >= H)[0])[0] != 0) and (np.where(D >= K)[0][0] >= np.where(mfp >= H)[0][0])):
        hom_id = 'homopause was found, but it is above the exobase'
    else:
        hom_id = np.where(D >= K)[0][0]
    
    return(hom_id)


def exobase_index_finder(H, mfp):
    """finds and returns the exobase index (where mfp becomes larger than H) if it exists"""

    if((np.shape(np.where(mfp >= H)[0])[0] == 0) or (np.shape(np.where(mfp >= H)[0])[0] == 1)):
        exo_id = 'exobase not found'
    else:
        exo_id = np.where(mfp >= H)[0][0]
    
    return(exo_id)


def Bates_extension(homopause_index, T_inf, M, R, T, p, z, rho, mmw, VMR, Kzz):
    """extends the atmosphere as a Bates profile up to 10^-16 bar given the value arrays, planet mass and radius, the homopause 
    index (from where the extenstion begins) and the heated up 'infinity' temperature at the very top of the atmosphere"""
    
    T = T[:homopause_index+1]
    p = p[:homopause_index+1]
    z = z[:homopause_index+1]
    rho = rho[:homopause_index+1]
    mmw = mmw[:homopause_index+1]
    VMR = {k: v.copy() for k, v in VMR.items()} #make a cope to not overwrite the dict
    for element in VMR:
        VMR[element] = VMR[element][:homopause_index+1]  
    Kzz = Kzz[:homopause_index+1]
        
    #extending p
    steps = 30
    p_ext = np.logspace(np.log10(p[-1]), -12, steps+1)[1:]  
        
    #extending the constant values: mmw, vmr, Kzz
    mmw_ext = np.full(steps, mmw[-1])
    mmw = np.concatenate((mmw, mmw_ext), axis=0)
    for element in VMR:
        last_VMR = VMR[element][-1]
        VMR_ext = np.full(steps, last_VMR)
        VMR[element] = np.concatenate((VMR[element], VMR_ext), axis=0)
    Kzz_ext = np.full(steps, Kzz[-1])
    Kzz = np.concatenate((Kzz, Kzz_ext), axis=0)
    
    #extending T     
    T_ext = T[-1] + ((T_inf - T[-1]) * (1-((p_ext/p[-1])**0.75)))
    T = np.concatenate((T, T_ext), axis=0)

    #extending density (ideal gas)
    rho_ext = ((mmw_ext/N_A)*(p_ext))/(k*T_ext)
    rho = np.concatenate((rho, rho_ext), axis=0)
    
    #extending height assuming constant g
    g_const = (G*M)/((R+z[-1])**2)
    H_arr = (k*T_ext)/((mmw_ext/N_A)*g_const)
    
    z_ext = z[-1] + (H_arr * np.log(p[-1]/p_ext))
    z = np.concatenate((z, z_ext), axis=0)
    
    p = np.concatenate((p, p_ext), axis=0)
    
    return(T, p, z, rho, mmw, VMR, Kzz)

    