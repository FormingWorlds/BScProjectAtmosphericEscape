import numpy as np
from scipy.constants import k, G, pi
import pandas as pd 

pd.set_option("display.float_format", "{:.3e}".format)

# core equations functions
def thermal_velocity(T, m):
    return np.sqrt(2 * k * T / m)

def jeans_parameter(M, m, T, r):
    return G * M * m / (k * T * r)

def effusion_velocity(v_th, lam):
    return (v_th / (2 * np.sqrt(pi))) * (1 + lam) * np.exp(-lam)


# average molecular mass for exobase calculation
def mean_mass(species, masses):
    n_stack = np.array(list(species.values()))
    m_array = np.array([masses[s] for s in species])
    
    return (n_stack * m_array[:, None]).sum(axis=0) / n_stack.sum(axis=0)
