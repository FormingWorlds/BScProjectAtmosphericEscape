import os
import numpy as np
import pandas as pd
from scipy.constants import k, atomic_mass

from jeans import jeans_escape
from constants import *
from extend_profile import *

def read_proteus_profile(path, R_planet, min_vmr=1e-21):

    df = pd.read_csv(path, sep="\t")

    # make sure profile goes bottom to top
    df = df.sort_values("Height [m]").reset_index(drop=True)

    r = R_planet + df["Height [m]"].to_numpy()
    T = df["Temperature [K]"].to_numpy()
    P = df["Pressure [Pa]"].to_numpy()

    n_tot = P / (k * T)
    #print(n_tot)

    species = {}

    for col in df.columns:
        if "[VMR]" not in col:
            continue

        sp = col.replace(" [VMR]", "")
        if sp not in SPECIES_MASSES:
            continue

        vmr = df[col].to_numpy()
        
        if np.nanmax(vmr) < min_vmr:
            continue

        species[sp] = vmr * n_tot
        #print(species[sp])

    return r, T, species, df


def read_bulk_properties(path, case):
    bulk = pd.read_csv(path, sep="\t")
    row = bulk[bulk["Case"] == case].iloc[0]
    return {
        "R_obs_m": row["R_obs [m]"],
        "R_int_m": row["R_int [m]"],
        "M_planet_kg": row["M_planet [kg]"],
        "F_xuv_W_m2": row["F_xuv [W/m2]"],
        "F_ins_W_m2": row["F_ins [W/m2]"],
        "MMW_g_mol": row["MMW [g/mol]"],
    }

def run_one_file(path, bulk):

    M_planet = bulk["M_planet_kg"]
    R_planet = bulk["R_int_m"]
    r, T, species, df = read_proteus_profile(path, R_planet)
    
    try:
        result = jeans_escape(
            r=r,
            T=T,
            species=species,
            species_masses=SPECIES_MASSES,
            M=M_planet,
            sigma="weighted",
            dayside=True,
        )
        
        profile_extended = False

    except ValueError:
        print("No exobase found, extending profile")

    r, T, species = extend_profile_exobase(
        r=r,
        T=T,
        species=species,
        species_masses = SPECIES_MASSES,
        M_planet = M_planet,
        z_extra=10e7,
        n_extra=10000,
    )
    
    result = jeans_escape(
            r=r,
            T=T,
            species=species,
            species_masses=SPECIES_MASSES,
            M=M_planet,
            sigma= "weighted",
            dayside=True,
        )

    profile_extended = True


    rows = []

    for sp, res in result["results"].items():
        rows.append({
            "file": os.path.basename(path),
            "weighted_mass_loss_kg_s": result["weighted_mass_loss_kg_s"],
            "species": sp,
            "Mdot_kg_s": res["Mdot (kg/s)"],
            "lambda_j": res["lambda_j"],
            "v_th_m_s": res["v_th (m/s)"],
            "effusion_velocity_m_s": res["effusion_velocity (m/s)"],
            "n_exo_m3": res["number_density (m^3)"],
            "exobase_altitude_km": result["exobase_altitude"] / 1e3,
            "exobase_radius_m": result["exobase_radius"],
            "T_exo_K": result["temperature"],
            "exobase_index": result["exobase_index"],
        })

    return rows

def summarize_case(rows):
    """
    Summarise all species escape rates for one atmosphere file.
    rows = list of dictionaries returned by run_one_file()
    """

    if len(rows) == 0:
        return None

    dominant = max(rows, key=lambda row: row["Mdot_kg_s"])

    summary = {
        "file": rows[0]["file"],
        "atmosphere_type": rows[0]["atmosphere_type"],
        "mass_case": rows[0]["mass_case"],
        "flux_case": rows[0]["flux_case"],

        "weighted_mass_loss_kg_s": rows[0]["weighted_mass_loss_kg_s"],
        "dominant_escaping_species": dominant["species"],
        "dominant_species_Mdot_kg_s": dominant["Mdot_kg_s"],

        "T_exo_K": rows[0]["T_exo_K"],
        "exobase_altitude_km": rows[0]["exobase_altitude_km"],
        "exobase_radius_m": rows[0]["exobase_radius_m"],
        "exobase_index": rows[0]["exobase_index"],
    }

    return summary

proteus_atmospheres = []

summary_rows = []

for comp in ['H2', 'H2O', 'CO2', 'N2']:

    for M in ['1_M_earth', '10_M_earth']:

        bulk_path = (
            f"PROTEUS data/"
            f"{comp}_atmospheres/"
            f"planet_bulk_properties_{comp}_atmospheres_{M}.csv"
        )

        try:

            for case_name in ["1_F_earth", "1000_F_earth"]:

                bulk = read_bulk_properties(
                    bulk_path,
                    case_name
                )

                profile_path = (
                    f"PROTEUS data/"
                    f"{comp}_atmospheres/"
                    f"{case_name}/"
                    f"{comp}_atmosphere_{M}_{case_name}.csv"
                )

                rows = run_one_file(
                    profile_path,
                    bulk,
                )
                # add metadata labels
                for row in rows:

                    row["atmosphere_type"] = comp
                    row["mass_case"] = M
                    row["flux_case"] = case_name

                proteus_atmospheres.extend(rows)
                summary_rows.append(summarize_case(rows))

                print(f"Finished: {comp} {M} {case_name}")

        except FileNotFoundError:

            print(f"Missing file for {comp} {M} {case_name}")

        except ValueError as e:

            print(f"Skipping {comp} {M}: {e} {case_name}")

results = pd.DataFrame(proteus_atmospheres)

results.to_csv(
    "proteus_jeans_escape_results.csv",
    index=False
)

print(results.head())

summary_results = pd.DataFrame(summary_rows)

summary_results.to_csv(
    "proteus_jeans_case_summary.csv",
    index=False
)
