import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.constants import k

from constants import SPECIES_MASSES
from proteus_escape import read_proteus_profile, read_bulk_properties
from extend_profile import extend_profile_exobase


OUTDIR = "Plots/PROTEUS_PT_extensions"
os.makedirs(OUTDIR, exist_ok=True)

atm_archetype = ["H2", "H2O", "CO2", "N2"]
masses = ["1_M_earth", "10_M_earth"]
fluxes = ["1_F_earth", "1000_F_earth"]

T_inf_values = [200, 300, 500, 1000, 2000, 3000, 5000, 7000]


def compute_pressure(T, species):
    """
    Compute pressure from total number density.

    P = n_tot k T
    """
    n_tot = np.sum(np.array(list(species.values())), axis=0)
    return n_tot * k * T


def plot_one_PT_extension(comp, mass_case, flux_case):
    bulk_path = (
        f"PROTEUS data/"
        f"{comp}_atmospheres/"
        f"planet_bulk_properties_{comp}_atmospheres_{mass_case}.csv"
    )

    profile_path = (
        f"PROTEUS data/"
        f"{comp}_atmospheres/"
        f"{flux_case}/"
        f"{comp}_atmosphere_{mass_case}_{flux_case}.csv"
    )

    if not os.path.exists(profile_path):
        print(f"Missing file: {profile_path}")
        return

    bulk = read_bulk_properties(bulk_path, flux_case)
    R_planet = bulk["R_int_m"]
    M_planet = bulk["M_planet_kg"]

    r, T, species, df = read_proteus_profile(profile_path, R_planet)

    P_original = compute_pressure(T, species)

    plt.figure(figsize=(5.5, 6.5))

    # Original PROTEUS profile
    plt.plot(
        T,
        P_original * 1e-5,
        color="black",
        linewidth=2.2,
        label="PROTEUS profile",
    )

    # Bates extensions
    for T_inf in T_inf_values:
        try:
            r_ext, T_ext, species_ext = extend_profile_exobase(
                r=r,
                T=T,
                species=species,
                species_masses=SPECIES_MASSES,
                M_planet=M_planet,
                T_inf=T_inf,
            )

            P_ext = compute_pressure(T_ext, species_ext)

            plt.plot(
                T_ext,
                P_ext * 1e-5,
                linewidth=1.4,
                label=rf"$T_\infty={T_inf}$ K",
            )

        except ValueError as e:
            print(f"Skipping {comp} {mass_case} {flux_case} T_inf={T_inf}: {e}")

    plt.yscale("log")
    plt.gca().invert_yaxis()

    plt.xlabel("Temperature [K]", fontsize=13)
    plt.ylabel("Pressure [bar]", fontsize=13)

    title = (
        f"{comp} atmosphere, "
        f"{mass_case.replace('_', ' ')}, "
        f"{flux_case.replace('_', ' ')}"
    )
    plt.title(title, fontsize=13)

    plt.legend(fontsize=8)
    plt.tight_layout()

    filename = f"PT_extension_{comp}_{mass_case}_{flux_case}.png"
    plt.savefig(os.path.join(OUTDIR, filename), dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved: {filename}")


for comp in atm_archetype:
    for mass_case in masses:
        for flux_case in fluxes:
            plot_one_PT_extension(comp, mass_case, flux_case)