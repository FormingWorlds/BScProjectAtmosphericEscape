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