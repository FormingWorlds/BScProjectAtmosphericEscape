import numpy as np
from constants import *


def effective_cross_section(species):
    n_tot = np.sum(list(species.values()), axis=0)
    sigma_eff = np.zeros_like(n_tot)

    for sp_i, n_i in species.items():
        if sp_i not in PARTICLE_RADII:
            continue

        r_i = PARTICLE_RADII[sp_i]
        X_i = n_i / n_tot

        for sp_j, n_j in species.items():
            if sp_j not in PARTICLE_RADII:
                continue

            X_j = n_j / n_tot

            r_j = PARTICLE_RADII[sp_j]

            sigma_ij = np.pi * (r_i + r_j)**2

            sigma_eff += X_i * X_j * sigma_ij

    return sigma_eff