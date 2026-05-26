import numpy as np
from scipy.constants import G, k
from scipy.integrate import cumulative_trapezoid

from exobase import find_exobase
from cross_section import effective_cross_section
from physics import mean_mass


def bates_extension(zeta, T0, T_inf, beta = 0.75):
    """extends using bates profile
    T(zeta) = T0 + (T_inf - T0) * (1 - exp(-beta*zeta))"""
    return T0 + (T_inf - T0) * (1.0 - np.exp(-beta * zeta))

def extend_profile_exobase(
    r,
    T,
    species,
    species_masses,
    M_planet,
    T_inf,
    zeta_extra=0.25,
    n_extra=10,
    zeta_max=30,
    ):
    """extend profile upwards to exobase using the Bates T profile
    Assumptions:
    - Bates temperature profile above PROTEUS top
    - pressure coordinate zeta = -ln(P/P_top)
    - constant VMRs from the top of the PROTEUS grid
    - hydrostatic relation dr/dzeta = H"""

    T0 = T[-1]
    r_current = r.copy()
    T_current = T.copy()
    species_current = {sp: n.copy() for sp, n in species.items()}

    zeta_added = 0

    while zeta_added < zeta_max:

        r_top = r_current[-1]
        T_top = T_current[-1]

        n_top = np.sum([n[-1] for n in species_current.values()])

        if n_top <= 0 or not np.isfinite(n_top):
            raise ValueError("Invalid top density before Bates extension.")

        X_top = {
            sp: species_current[sp][-1] / n_top
            for sp in species_current
        }

        m_mean_top = sum(
            X_top[sp] * species_masses[sp]
            for sp in species_current
        )

        zeta_grid = np.linspace(0, zeta_extra, n_extra)

        zeta_total = zeta_grid + zeta_added

        T_ext = bates_extension(zeta_total, T_top, T_inf, beta=0.75)
        
        g_ext = G * M_planet / r_top**2
        H_ext = k * T_ext / (m_mean_top * g_ext)

        dr_ext = cumulative_trapezoid(
            H_ext,
            zeta_grid,
            initial=0,
        )

        r_ext =  r_top + dr_ext

        if r_ext[-1] > 1e12:   # more than ~7000 R_earth, unphysical?
            raise ValueError("Extension reached unphysical radius — exobase not found")
        
        n_extended = n_top * np.exp(-zeta_grid) * (T_top/T_ext)

        r_ext = r_ext[1:]    # remove duplicate r_top
        n_ext_tot = n_extended[1:]     # remove corresponding density
        T_ext = T_ext[1:]

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
            idx = find_exobase(
                r_new,
                n_tot_new,
                T_new,
                m_mean_new,
                M_planet,
                sigma_new,
            )

            r_final = r_new[:idx + 1]
            T_final = T_new[:idx + 1]

            species_final = {
                sp: species_new[sp][:idx + 1]
                for sp in species_new
            }

            return r_final, T_final, species_final

        except ValueError:
            r_current = r_new
            T_current = T_new
            species_current = species_new
            zeta_added += zeta_extra

    raise ValueError("No exobase found after extending profile with Bates")