import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

outdir = "Plots/Extension_diagnostics"
os.makedirs(outdir, exist_ok=True)

summary = pd.read_csv("proteus_jeans_case_summary.csv")

atm_archetype = ["CO2", "N2", "H2O", "H2"]
mass_cases = ["1_M_earth", "10_M_earth"]
flux_cases = ["1_F_earth", "1000_F_earth"]

species_colors = {
    "CO2": "#E69F00",
    "H2":  "#332288",
    "H2O": "#56B4E9",
    "N2":  "#882255",
}

mass_labels = {
    "1_M_earth": r"$1\,M_{\oplus}$",
    "10_M_earth": r"$10\,M_{\oplus}$",
}

flux_labels = {
    "1_F_earth": r"1 $F_{\oplus}$",
    "1000_F_earth": r"1000 $F_{\oplus}$",
}

flux_markers = {
    "1_F_earth": "o",
    "1000_F_earth": "^",
}

# If your exobase_altitude_km is measured from the profile bottom,
# and the PROTEUS top is the last height of the original file,
# then extension height is:
# extension_height = exobase altitude - PROTEUS top height.
#
# To compute it, we need to read the original profiles again.

def get_profile_top_height_km(row):
    comp = row["atmosphere_type"]
    mass = row["mass_case"]
    flux = row["flux_case"]

    path = (
        f"PROTEUS data/"
        f"{comp}_atmospheres/"
        f"{flux}/"
        f"{comp}_atmosphere_{mass}_{flux}.csv"
    )

    df = pd.read_csv(path, sep="\t")
    return df["Height [m]"].max() / 1e3


summary["profile_top_height_km"] = summary.apply(get_profile_top_height_km, axis=1)

summary["extension_height_km"] = (
    summary["exobase_altitude_km"] -
    summary["profile_top_height_km"]
)

summary.to_csv("extension_diagnostics_summary.csv", index=False)

atm_handles = [
    Line2D([0], [0], color=species_colors[a], lw=2, label=a)
    for a in atm_archetype
]

flux_handles = [
    Line2D([0], [0], marker=flux_markers[f], color="black",
           linestyle="None", label=flux_labels[f])
    for f in flux_cases
]

invalid_handle = Line2D(
    [0], [0],
    marker="x",
    color="black",
    linestyle="None",
    markersize=8,
    markeredgewidth=2,
    label=r"Hydrodynamic onset ($\lambda_J<1.5$)",
)

all_handles = atm_handles + flux_handles + [invalid_handle]



fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for idx, mass in enumerate(mass_cases):
    ax = axes[idx]

    for atm in atm_archetype:
        for flux in flux_cases:

            s = summary[
                (summary["mass_case"] == mass) &
                (summary["atmosphere_type"] == atm) &
                (summary["flux_case"] == flux)
            ].sort_values("T_inf")

            if s.empty:
                continue

            valid = s[s["jeans_valid"]]
            invalid = s[~s["jeans_valid"]]

            ax.plot(
                valid["T_inf"],
                valid["extension_height_km"],
                color=species_colors[atm],
                linewidth=1.5,
                alpha=0.75,
            )

            ax.scatter(
                valid["T_inf"],
                valid["extension_height_km"],
                color=species_colors[atm],
                marker=flux_markers[flux],
                s=55,
            )

            ax.scatter(
                invalid["T_inf"],
                invalid["extension_height_km"],
                color=species_colors[atm],
                marker="x",
                s=80,
                linewidths=2,
            )

    ax.axhline(1, color="grey", linestyle="--", linewidth=1)
    ax.set_yscale("log")

    ax.set_xlabel(r"$T_{\infty}$ [K]", fontsize=13)
    ax.set_title(mass_labels[mass], fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)

axes[0].set_ylabel("Extension height to exobase [km]", fontsize=13)

fig.legend(
    handles=all_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.12),
    ncol=4,
    frameon=False,
    fontsize=9,
)

plt.tight_layout(rect=[0, 0.12, 1, 1])
plt.savefig(os.path.join(outdir, "extension_height_vs_Tinf.png"),
            dpi=300, bbox_inches="tight")
plt.close()



from scipy.constants import k, G
from scipy.integrate import cumulative_trapezoid

from extend_profile import bates_extension
from constants import SPECIES_MASSES

outdir = "Plots/Bates_extension"
os.makedirs(outdir, exist_ok=True)

T_inf_values = [200, 300, 500, 1000, 2000, 3000, 5000, 7000]

atm_archetype = ["CO2", "N2", "H2O", "H2"]
mass_case = "10_M_earth"
flux_case = "1000_F_earth"

def read_bulk(comp, mass_case, flux_case):
    path = (
        f"PROTEUS data/{comp}_atmospheres/"
        f"planet_bulk_properties_{comp}_atmospheres_{mass_case}.csv"
    )
    bulk = pd.read_csv(path, sep="\t")
    row = bulk[bulk["Case"] == flux_case].iloc[0]
    return row["R_int [m]"], row["M_planet [kg]"]

def plot_bates_T_vs_altitude(comp):
    profile_path = (
        f"PROTEUS data/{comp}_atmospheres/{flux_case}/"
        f"{comp}_atmosphere_{mass_case}_{flux_case}.csv"
    )

    R_planet, M_planet = read_bulk(comp, mass_case, flux_case)

    df = pd.read_csv(profile_path, sep="\t")
    df = df.sort_values("Height [m]").reset_index(drop=True)

    height = df["Height [m]"].to_numpy()
    T = df["Temperature [K]"].to_numpy()
    P = df["Pressure [Pa]"].to_numpy()

    r = R_planet + height
    n_tot = P / (k * T)

    species = {}
    for col in df.columns:
        if "[VMR]" not in col:
            continue

        sp = col.replace(" [VMR]", "")

        if sp not in SPECIES_MASSES:
            continue

        vmr = df[col].to_numpy()

        if np.nanmax(vmr) < 1e-21:
            continue

        species[sp] = vmr * n_tot

    # top of PROTEUS profile
    T0 = T[-1]
    r_top = r[-1]
    height_top = height[-1]
    n_top = sum(n[-1] for n in species.values())

    X_top = {
        sp: species[sp][-1] / n_top
        for sp in species
    }

    m_mean_top = sum(
        X_top[sp] * SPECIES_MASSES[sp]
        for sp in species
    )

    P_top = n_top * k * T0
    P_target = 1e-12

    zeta_max = np.log(P_top / P_target)
    zeta = np.linspace(0, zeta_max, 1000)

    g_top = G * M_planet / r_top**2

    plt.figure(figsize=(7, 5))

    # original PROTEUS profile
    plt.plot(
        T,
        height / 1e3,
        color="black",
        linewidth=2.5,
        label="PROTEUS profile"
    )

    for T_inf in T_inf_values:
        T_ext = bates_extension(zeta, T0, T_inf, beta=0.75)

        H_ext = k * T_ext / (m_mean_top * g_top)

        dz_ext = cumulative_trapezoid(
            H_ext,
            zeta,
            initial=0
        )

        height_ext = height_top + dz_ext

        plt.plot(
            T_ext,
            height_ext / 1e3,
            linewidth=1.7,
            label=fr"$T_\infty={T_inf}$ K"
        )

    plt.xlabel("Temperature [K]")
    plt.ylabel("Altitude [km]")
    plt.title(
        fr"{comp} atmosphere, "
        fr"$10\,M_\oplus$, $1000\,F_\oplus$"
    )

    plt.yscale("log")
    plt.legend(fontsize=8)
    plt.tight_layout()

    plt.savefig(
        os.path.join(outdir, f"{comp}_10M_1000F_Bates_T_vs_altitude.png"),
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()


for comp in atm_archetype:
    plot_bates_T_vs_altitude(comp)


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.constants import k

from extend_profile import bates_extension
from constants import SPECIES_MASSES

outdir = "Plots/Bates_extension"
os.makedirs(outdir, exist_ok=True)

T_inf_values = [200, 300, 500, 1000, 2000, 3000, 5000, 7000]

comp = "H2O"              # choose one representative atmosphere
mass_case = "10_M_earth"
flux_case = "1000_F_earth"

profile_path = (
    f"PROTEUS data/{comp}_atmospheres/{flux_case}/"
    f"{comp}_atmosphere_{mass_case}_{flux_case}.csv"
)

df = pd.read_csv(profile_path, sep="\t")
df = df.sort_values("Height [m]").reset_index(drop=True)

T0 = df["Temperature [K]"].iloc[-1]
P_top = df["Pressure [Pa]"].iloc[-1]

P_target = 1e-12  # Pa

P_ext = np.logspace(
    np.log10(P_top),
    np.log10(P_target),
    1000
)

zeta = -np.log(P_ext / P_top)

plt.figure(figsize=(7, 5))

for T_inf in T_inf_values:
    T_ext = bates_extension(
        zeta,
        T0,
        T_inf,
        beta=0.75
    )

    plt.plot(
        T_ext,
        zeta,
        linewidth=1.7,
        label=fr"$T_\infty={T_inf}$ K"
    )

plt.xlabel("Temperature [K]")
plt.ylabel(r"$\zeta = -\ln(P/P_{\rm top})$")
plt.title(
    fr"{comp} atmosphere, "
    fr"$10\,M_\oplus$, $1000\,F_\oplus$"
)

plt.legend(fontsize=8, loc = "lower right")
plt.tight_layout()

plt.savefig(
    os.path.join(outdir, f"{comp}_10M_1000F_Bates_T_vs_zeta_pressure.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.close()