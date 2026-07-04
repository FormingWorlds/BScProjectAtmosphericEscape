import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

outdir = "Plots/Dissociation_comparison"
os.makedirs(outdir, exist_ok=True)

T_choice = 5000
mass_choice = "10_M_earth"
flux_choice = "1000_F_earth"

atm_order = ["CO2", "N2", "H2O", "H2"]

atm_labels = {
    "CO2": r"CO$_2$",
    "N2": r"N$_2$",
    "H2O": r"H$_2$O",
    "H2": r"H$_2$",
}

constant = pd.read_csv("Outputs/proteus_jeans_case_summary.csv")
diss = pd.read_csv("Outputs/proteus_jeans_case_summary_dissociation.csv")

for df in [constant, diss]:
    if "unphysical_extension" not in df.columns:
        df["unphysical_extension"] = False

constant["case"] = "Constant VMR"
diss["case"] = "Full dissociation"

df = pd.concat([constant, diss], ignore_index=True)

s = df[
    (df["T_inf"] == T_choice) &
    (df["mass_case"] == mass_choice) &
    (df["flux_case"] == flux_choice) &
    (~df["unphysical_extension"])
].copy()

s["atmosphere_type"] = pd.Categorical(
    s["atmosphere_type"],
    categories=atm_order,
    ordered=True
)

s = s.sort_values(["atmosphere_type", "case"])


# grouped bar chart

pivot = s.pivot(
    index="atmosphere_type",
    columns="case",
    values="dominant_species_Mdot_kg_s"
).reindex(atm_order)

x = np.arange(len(atm_order))
width = 0.36

fig, ax = plt.subplots(figsize=(7, 5))

ax.bar(
    x - width / 2,
    np.log10(pivot["Constant VMR"]),
    width,
    label="Constant VMR",
    color="#882255"
)

ax.bar(
    x + width / 2,
    np.log10(pivot["Full dissociation"]),
    width,
    label="Full dissociation",
    color="#56B4E9",
)

for i, atm in enumerate(atm_order):

    for case_name, x_offset in [
        ("Constant VMR", -width / 2),
        ("Full dissociation", width / 2),
    ]:
        selected = s[
            (s["atmosphere_type"] == atm) &
            (s["case"] == case_name)
        ]

        if selected.empty:
            continue

        val = selected["dominant_species_Mdot_kg_s"].iloc[0]
        sp = selected["dominant_escaping_species"].iloc[0]

        if pd.isna(val) or val <= 0:
            continue

        ax.text(
            x[i] + x_offset,
            np.log10(val) + 0.05,
            sp,
            ha="center",
            va="bottom",
            fontsize=9,
        )

ax.set_xticks(x)
ax.set_xticklabels([atm_labels[a] for a in atm_order])

ax.set_ylabel(
    r"$\log_{10}(\dot{M}_{\rm dominant})$ [kg s$^{-1}$]"
)
ax.set_xlabel("Atmosphere type")

ax.set_title(
    fr"Dominant escape rate at $T_\infty={T_choice}$ K, "
    fr"$10\,M_\oplus$, $1000\,F_\oplus$"
)

ax.legend(frameon=False)
plt.tight_layout()

plt.savefig(
    os.path.join(outdir, "dominant_Mdot_bar_constant_vs_dissociation.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.close()


cols = [
    "atmosphere_type", "mass_case", "flux_case", "T_inf",
    "species", "Mdot_kg_s", "lambda_j", "n_exo_m3",
    "exobase_altitude_km", "exobase_radius_m", "T_exo_K"
]

df = pd.read_csv("Outputs/proteus_jeans_escape_results_dissociation.csv")

case = df[
    (df["atmosphere_type"] == "N2") &
    (df["flux_case"] == "1000_F_earth") &
    (df["T_inf"] == 5000) &
    (~df["unphysical_extension"])
]

print(
    case[cols]
    .sort_values(["mass_case", "Mdot_kg_s"], ascending=[True, False])
    .to_string(index=False)
)