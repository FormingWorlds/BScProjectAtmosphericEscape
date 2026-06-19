import numpy as np
import pandas as pd
from scipy.constants import k
from scipy.integrate import cumulative_trapezoid
from scipy.constants import G

from constants import SPECIES_MASSES, DISSOCIATION
from extend_dissociation import bates_extension, apply_full_photodissociation

# Choose representative case
comp = "H2O"
mass_case = "10_M_earth"
flux_case = "1000_F_earth"
T_inf = 1000
P_target = 1e-12
min_vmr = 1e-21

profile_path = (
    f"PROTEUS data/{comp}_atmospheres/{flux_case}/"
    f"{comp}_atmosphere_{mass_case}_{flux_case}.csv"
)

bulk_path = (
    f"PROTEUS data/{comp}_atmospheres/"
    f"planet_bulk_properties_{comp}_atmospheres_{mass_case}.csv"
)

bulk = pd.read_csv(bulk_path, sep="\t")
bulk_row = bulk[bulk["Case"] == flux_case].iloc[0]

R_planet = bulk_row["R_int [m]"]
M_planet = bulk_row["M_planet [kg]"]

df = pd.read_csv(profile_path, sep="\t")
df = df.sort_values("Height [m]").reset_index(drop=True)

r = R_planet + df["Height [m]"].to_numpy()
T = df["Temperature [K]"].to_numpy()
P = df["Pressure [Pa]"].to_numpy()

n_tot = P / (k * T)

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

# Top of PROTEUS profile
T0 = T[-1]
r_top = r[-1]
T_top = T[-1]
n_top = sum(n[-1] for n in species.values())

X_top = {
    sp: species[sp][-1] / n_top
    for sp in species
}

m_mean_top = sum(
    X_top[sp] * SPECIES_MASSES[sp]
    for sp in species
)

P_top = n_top * k * T_top

zeta_max = np.log(P_top / P_target)
zeta = np.linspace(0, zeta_max, 1000)

T_ext = bates_extension(zeta, T0, T_inf, beta=0.75)

P_ext = P_top * np.exp(-zeta)
n_ext_tot = P_ext / (k * T_ext)

species_ext_const = {
    sp: X_top[sp] * n_ext_tot
    for sp in species
}

species_ext_diss = apply_full_photodissociation(species_ext_const)

# Diagnostic comparison at top of extension and upper end of extension
rows = []

all_species = sorted(
    set(species_ext_const.keys()) | set(species_ext_diss.keys())
)

for sp in all_species:
    n_const_start = species_ext_const.get(sp, np.zeros_like(zeta))[0]
    n_diss_start = species_ext_diss.get(sp, np.zeros_like(zeta))[0]

    n_const_end = species_ext_const.get(sp, np.zeros_like(zeta))[-1]
    n_diss_end = species_ext_diss.get(sp, np.zeros_like(zeta))[-1]

    rows.append({
        "species": sp,
        "n_const_start_m3": n_const_start,
        "n_diss_start_m3": n_diss_start,
        "ratio_start_diss_const": (
            n_diss_start / n_const_start
            if n_const_start > 0 else np.nan
        ),
        "n_const_end_m3": n_const_end,
        "n_diss_end_m3": n_diss_end,
        "ratio_end_diss_const": (
            n_diss_end / n_const_end
            if n_const_end > 0 else np.nan
        ),
    })

diagnostic = pd.DataFrame(rows)

diagnostic.to_csv(
    f"dissociation_diagnostic_{comp}_{mass_case}_{flux_case}_Tinf{T_inf}.csv",
    index=False
)

pd.set_option("display.max_rows", None)
pd.set_option("display.float_format", "{:.3e}".format)

print()
print(f"Case: {comp}, {mass_case}, {flux_case}, T_inf={T_inf} K")
print(f"P_top = {P_top:.3e} Pa")
print(f"T_top = {T_top:.3e} K")
print(f"n_top = {n_top:.3e} m^-3")
print()

print(diagnostic)

# Print especially relevant species
print()
print("Selected species:")
for sp in ["H", "H2", "H2O", "O", "OH", "CO2", "N2"]:
    if sp in all_species:
        row = diagnostic[diagnostic["species"] == sp].iloc[0]
        print(
            f"{sp:4s} | "
            f"const start = {row['n_const_start_m3']:.3e}, "
            f"diss start = {row['n_diss_start_m3']:.3e}, "
            f"ratio = {row['ratio_start_diss_const']:.3e}"
        )