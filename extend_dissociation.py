import numpy as np
from scipy.constants import G, k
from scipy.integrate import cumulative_trapezoid

from exobase import find_exobase
from cross_section import effective_cross_section
from physics import mean_mass
from constants import DISSOCIATION


def apply_full_photodissociation(species_ext):
    """
    Fully dissociate selected molecules in the extended region.
    """

    # Make a copy first to not ruin initial data
    species_out = {
        sp: n.copy()
        for sp, n in species_ext.items()
    }

    for parent, products in DISSOCIATION.items():

        # If this molecule is not present, skip it
        if parent not in species_out:
            continue

        # keep initial density stored
        n_parent = species_out[parent].copy()

        # Fully remove parent molecule from extension because it dissociated
        species_out[parent] = np.zeros_like(n_parent)

        # Add products
        for product, value in products.items():

            # If product does not exist yet, create empty array
            if product not in species_out:
                species_out[product] = np.zeros_like(n_parent)

            # Add the amount of product
            species_out[product] += value * n_parent

    return species_out


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
    ):
    """extend profile upwards to exobase using the Bates T profile,
    upper limit full dissociation case
    Assumptions:
    - Bates temperature profile above PROTEUS top
    - pressure coordinate zeta = -ln(P/P_top)
    - full dissociation in extension
    - hydrostatic relation dr/dzeta = H"""

    T0 = T[-1]
    r_current = r.copy()
    T_current = T.copy()
    species_current = {sp: n.copy() for sp, n in species.items()}

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

    P_target = 1e-12          # bar, for space like pressure
    P_top = n_top * k * T_top

    zeta_max = np.log(P_top / P_target)
    zeta = np.linspace(0, zeta_max, 1000)

    T_ext = bates_extension(zeta, T0, T_inf, beta=0.75)

    P_ext= P_top * np.exp(-zeta)

    n_ext_tot = P_ext / (k * T_ext)

    g_ext = G * M_planet / r_top**2
    H_ext = k * T_ext / (m_mean_top * g_ext)

    dr_ext = cumulative_trapezoid(
        H_ext,
        zeta,
        initial=0,
    )

    r_ext =  r_top + dr_ext

    if r_ext[-1] > 20 * r[0]:   # more than ~20 Rp, unphysical
        raise ValueError("UNPHYSICAL EXTENSION")

    n_extended = n_top * np.exp(-zeta) * (T_top/T_ext)

    r_ext = r_ext[1:]    # remove duplicate r_top
    n_ext_tot = n_extended[1:]     # remove corresponding density
    T_ext = T_ext[1:]

    species_ext = {
        sp: X_top[sp] * n_ext_tot
        for sp in species_current
    }

    species_ext = apply_full_photodissociation(species_ext) #change to dissociation instead of ct VMR

    r_new = np.concatenate([r_current, r_ext])
    T_new = np.concatenate([T_current, T_ext])

    #here we include species that might only appear after dissociation with 0 before extension

    all_species = set(species.keys()) | set(species_ext.keys()) #unite all species 

    species_new = {}

    for sp in all_species:
        original = species.get(sp, np.zeros_like(r))
        extension = species_ext.get(sp, np.zeros_like(r_ext))

        species_new[sp] = np.concatenate([original, extension])

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

    raise ValueError("No exobase found after extending profile with Bates")