import numpy as np
from scipy.constants import G, k
from scipy.integrate import cumulative_trapezoid

def extend_profile_exobase(
        r,
        T,
        species,
        species_masses,
        M_planet,
        z_extra=10e7,
        n_extra=10e5,
):
    """extend profile upwards to exobase assuming 
    hydrostatic equilibrium and isothermal extension, 
    constant VMR from top of grid"""
    r_top = r[-1]
    T_top = T[-1]

    n_top = np.sum([n[-1] for n in species.values()])

    r_ext = np.linspace(r_top, r_top +z_extra, n_extra)
    T_ext = np.full_like(r_ext, T_top)

    m_mean_top = sum(species[sp][-1] * species_masses[sp] for sp in species) / n_tot

    g_ext = G * M_planet / r_ext**2
    H_ext = k * T_ext / (m_mean_top * g_ext)

    ln_n = cumulative_trapezoid(
        -1 / H_ext,
        r_ext,
        intial = 0,
    )
    n_ext_tot = n_top * np.exp(ln_n)

    species_ext = {}

    for sp, n in species.items():
        species_ext[sp] = n[-1]/n_top *n_ext_tot
    
    r_new = np.concatenate ([r, r_ext])
    T_new = np.concatenate ([T, T_ext])


    species_new = {}

    for sp in species:
        species_new[sp] = np.concatenate([species[sp], species_ext[sp]])

    return r_new, T_new, species_new




