import numpy as np
from scipy.constants import k, G, N_A


def scale_height(T, R, z, M, mmw):
    """finds density scale height"""
    
    H = (k * T * ((R + z)**2))/(G * M * (mmw/N_A))
    return(H)