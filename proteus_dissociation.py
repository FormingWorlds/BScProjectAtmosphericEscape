import os
import numpy as np
import pandas as pd
from scipy.constants import k, atomic_mass

from jeans import jeans_escape
from constants import *
from extend_dissociation import *

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
    profile_extended = False

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
        

    except ValueError:
        profile_extended = True
        print("No exobase found, extending profile")

    if profile_extended == True: 
        T_inf_values = [200, 300, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000] #sensitivity testing
        rows = []
        for T_inf in T_inf_values:
            try:
                r_ext, T_ext, species_ext = extend_profile_exobase(
                    r=r,
                    T=T,
                    species=species,
                    species_masses = SPECIES_MASSES,
                    M_planet=M_planet,
                    T_inf=T_inf,
                    )
                
                result = jeans_escape(
                        r=r_ext,
                        T=T_ext,
                        species=species_ext,
                        species_masses=SPECIES_MASSES,
                        M=M_planet,
                        sigma= "weighted",
                        dayside=True,
                        )
            except ValueError as e:

                if str(e) == "UNPHYSICAL EXTENSION":

                    rows.append({
                        "file": os.path.basename(path),
                        "T_inf": T_inf,
                        "weighted_mass_loss_kg_s": np.nan,
                        "species": "None",
                        "Mdot_kg_s": np.nan,
                        "lambda_j": np.nan,
                        "v_th_m_s": np.nan,
                        "effusion_velocity_m_s": np.nan,
                        "n_exo_m3": np.nan,
                        "exobase_altitude_km": np.nan,
                        "exobase_radius_m": np.nan,
                        "T_exo_K": np.nan,
                        "exobase_index": np.nan,
                        "jeans_valid": False,
                        "unphysical_extension": True,
                        "escape_regime": "unphysical_extension",
                    })

                    print(
                        f"{os.path.basename(path)} "
                        f"for T_inf={T_inf}: {e}"
                    )
                    

                    continue

                else:
                    print(
                        f"Skipping {os.path.basename(path)} "
                        f"for T_inf={T_inf}: {e}"
                    )
                    continue

            for sp, res in result["results"].items():
                rows.append({
                    "file": os.path.basename(path),
                    "T_inf": T_inf,
                    "weighted_mass_loss_kg_s": result["weighted_mass_loss_kg_s"],
                    "species": sp,
                    "M_total_kg_s": result["total_mass_loss_kg_s"],
                    "Mdot_kg_s": res["Mdot (kg/s)"],
                    "lambda_j": res["lambda_j"],
                    "v_th_m_s": res["v_th (m/s)"],
                    "effusion_velocity_m_s": res["effusion_velocity (m/s)"],
                    "n_exo_m3": res["number_density (m^3)"],
                    "exobase_altitude_km": result["exobase_altitude"] / 1e3,
                    "exobase_radius_m": result["exobase_radius"],
                    "T_exo_K": result["exobase_temperature"],
                    "exobase_index": result["exobase_index"],
                    "jeans_valid": res["lambda_j"] >= 1.5,
                    "escape_regime": "Jeans" if res["lambda_j"] >= 1.5 else "hydrodynamic_candidate",
                    "unphysical_extension": False,
                    })

    return rows

def summarize_case(rows):
    """
    Summarise all species escape rates for one atmosphere file.
    rows = list of dictionaries returned by run_one_file()
    """
    valid_rows = [
    r for r in rows
    if not r["unphysical_extension"]
]

    if len(valid_rows) == 0:
        return {
            "file": rows[0]["file"],
            "T_inf": rows[0]["T_inf"],
            "atmosphere_type": rows[0]["atmosphere_type"],
            "mass_case": rows[0]["mass_case"],
            "flux_case": rows[0]["flux_case"],
            "unphysical_extension": True,
            "dominant_escaping_species": None,
            "dominant_species_Mdot_kg_s": np.nan,
            "dominant_lambda_j": np.nan,
            "weighted_mass_loss_kg_s": np.nan,
            "jeans_valid": False,
        }

    dominant = max(
        valid_rows,
        key=lambda r: r["Mdot_kg_s"]
    )
    if len(rows) == 0:
        return None

    summary = {
        "file": rows[0]["file"],
        "T_inf": rows[0]["T_inf"],
        "atmosphere_type": rows[0]["atmosphere_type"],
        "mass_case": rows[0]["mass_case"],
        "flux_case": rows[0]["flux_case"],

        "weighted_mass_loss_kg_s": rows[0]["weighted_mass_loss_kg_s"], #same for all species
        "M_total_kg_s": rows[0]["M_total_kg_s"],
        "dominant_escaping_species": dominant["species"],
        "dominant_lambda_j": dominant["lambda_j"],
        "dominant_species_Mdot_kg_s": dominant["Mdot_kg_s"],
        "jeans_valid": dominant["lambda_j"] >= 1.5,

        "T_exo_K": rows[0]["T_exo_K"],
        "exobase_altitude_km": rows[0]["exobase_altitude_km"],
        "exobase_radius_m": rows[0]["exobase_radius_m"],
        "exobase_index": rows[0]["exobase_index"],
        "unphysical_extension": rows[0]["unphysical_extension"],
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

                if len(rows) == 0:
                    print(f"No valid results for {comp} {M} {case_name}, skipping summary")
                    continue
                
                for T_inf, group in pd.DataFrame(rows).groupby("T_inf"):
                    summary_rows.append(summarize_case(group.to_dict("records")))

                print(f"Finished: {comp} {M} {case_name}")

        except FileNotFoundError:

            print(f"Missing file for {comp} {M} {case_name}")

        except ValueError as e:

            print(f"Skipping {comp} {M}: {e} {case_name}")

results = pd.DataFrame(proteus_atmospheres)

results.to_csv(
    "Outputs/proteus_jeans_escape_results_dissociation.csv",
    index=False
)

print(results.head())

summary_results = pd.DataFrame(summary_rows)

summary_results.to_csv(
    "Outputs/proteus_jeans_case_summary_dissociation.csv",
    index=False
)
