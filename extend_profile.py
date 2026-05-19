import numpy as np
from scipy.constants import G, k
from scipy.integrate import cumulative_trapezoid

from exobase import find_exobase
from cross_section import effective_cross_section
from physics import mean_mass


def extend_profile_exobase(
    r,
    T,
    species,
    species_masses,
    M_planet,
    z_extra=1e6,
    n_extra=10000,
    max_extra=1e9,
):
    """extend profile upwards to exobase assuming hydrostatic equilibrium and isothermal extension, 
    constant VMR from top of grid"""

    r_current = r.copy()
    T_current = T.copy()
    species_current = {sp: n.copy() for sp, n in species.items()}

    total_added = 0.0

    while total_added < max_extra:

        r_top = r_current[-1]
        T_top = T_current[-1]

        n_top = np.sum([n[-1] for n in species_current.values()])

        X_top = {
            sp: species_current[sp][-1] / n_top
            for sp in species_current
        }

        m_mean_top = sum(
            X_top[sp] * species_masses[sp]
            for sp in species_current
        )

        r_ext = np.linspace(r_top, r_top + z_extra, n_extra)[1:]
        T_ext = np.full_like(r_ext, T_top)

        g_ext = G * M_planet / r_ext**2
        H_ext = k * T_ext / (m_mean_top * g_ext)

        ln_n = cumulative_trapezoid(
            -1 / H_ext,
            r_ext,
            initial=0,
        )

        n_ext_tot = n_top * np.exp(ln_n)

        species_ext = {
            sp: X_top[sp] * n_ext_tot
            for sp in species_current
        }

        r_new = np.concatenate([r_current, r_ext])
        T_new = np.concatenate([T_current, T_ext])

        species_new = {
            sp: np.concatenate([species_current[sp], species_ext[sp]])
            for sp in species_current
        }

        n_tot_new = np.sum(
            np.array(list(species_new.values())),
            axis=0
        )

        m_mean_new = mean_mass(species_new, species_masses)
        sigma_new = effective_cross_section(species_new)

        try:
            find_exobase(
                r_new,
                n_tot_new,
                T_new,
                m_mean_new,
                M_planet,
                sigma_new,
            )

            return r_new, T_new, species_new

        except ValueError:
            r_current = r_new
            T_current = T_new
            species_current = species_new
            total_added += z_extra

    raise ValueError("No exobase found after extending profile")