import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

outdir = "Plots/Extension_diagnostics"
os.makedirs(outdir, exist_ok=True)

summary = pd.read_csv("Outputs/proteus_jeans_case_summary.csv")

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

# extension height is:
# extension_height = exobase altitude - PROTEUS top height.


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

summary.to_csv("Outputs/extension_diagnostics_summary.csv", index=False)

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



fig, axes = plt.subplots(2, 1, figsize=(5, 7), sharey=True)

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

            valid = s[s["jeans_valid"] & (~s["unphysical_extension"])]
            invalid = s[~s["jeans_valid"]& (~s["unphysical_extension"])]

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

    ax.set_xlabel(r"$T_{\infty}$ [K]", fontsize=11)
    ax.set_title(mass_labels[mass], fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)

axes[0].set_ylabel("Extension height to exobase [km]", fontsize=11)
axes[1].set_ylabel("Extension height to exobase [km]", fontsize=11)

fig.legend(
    handles=all_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.04),
    ncol=4,
    frameon=False,
    fontsize=9,
)

plt.tight_layout()
plt.savefig(os.path.join(outdir, "extension_height_vs_Tinf.png"),
            dpi=300, bbox_inches="tight")
plt.close()
