from scipy.constants import G, k
import numpy as np

# walking up to exobase
def find_exobase(r, n_tot, T, m_mean, M, sigma):
    g = G * M / r**2
    H = k * T / (m_mean * g)
    mfp = 1 / (sigma * n_tot)
    idx = np.where(mfp > H)[0]  
    if len(idx) == 0:
        raise ValueError("No exobase found")
    
    return idx[0]