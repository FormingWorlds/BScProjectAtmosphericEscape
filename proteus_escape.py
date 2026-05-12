import os
import glob
import numpy as np
import pandas as pd
from scipy.constants import k, atomic_mass

from jeans import jeans_escape


R_EARTH = 6.371e6
M_EARTH = 5.972e24

SPECIES_MASSES = {
    "H": 1 * atomic_mass,
    "H2": 2 * atomic_mass,
    "He": 4 * atomic_mass,
    "C": 12 * atomic_mass,
    "N": 14 * atomic_mass,
    "O": 16 * atomic_mass,
    "N2": 28 * atomic_mass,
    "O2": 32 * atomic_mass,
    "CO": 28 * atomic_mass,
    "CO2": 44 * atomic_mass,
    "H2O": 18 * atomic_mass,
    "OH": 17 * atomic_mass,
    "NO": 30 * atomic_mass,
    "NO2": 46 * atomic_mass,
    "N2O": 44 * atomic_mass,
    "CH4": 16 * atomic_mass,
    "NH3": 17 * atomic_mass,
}

def read_proteus_profile(path, R_planet):
    df = pd.read_csv(path, sep="\t")

    r = R_planet + df["Height [m]"].to_numpy()
    T = df["Temperature [K]"].to_numpy()
    P = df["Pressure [Pa]"].to_numpy()

    n_tot = P / (k * T)

    species = {}

    for col in df.columns:
        if not col.endswith("[VMR]"):
            continue

        sp = col.replace(" [VMR]", "")
        if sp not in SPECIES_MASSES:
            continue

        vmr = df[col].to_numpy()

        species[sp] = vmr * n_tot

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
    R_planet = bulk["R_obs_m"]
    r, T, species, df = read_proteus_profile(path, R_planet)

    result = jeans_escape(
        r=r,
        T=T,
        species=species,
        species_masses=SPECIES_MASSES,
        M=M_planet,
        sigma=sigma,
        dayside=True,
    )

    rows = []

    for sp, res in result["results"].items():
        rows.append({
            "file": os.path.basename(path),
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


