import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

os.makedirs("Plots/PROTEUS plots", exist_ok=True)

df = pd.read_csv("proteus_jeans_escape_results.csv")

case_cols = ["file", "T_inf", "atmosphere_type", "mass_case", "flux_case"]

atm_archetype = ["CO2", "N2", "H2O", "H2"]
instellation = ["1_F_earth", "1000_F_earth"]
mass = ["1_M_earth", "10_M_earth"]

species_colors = {
    # Atmosphere archetypes (distinctive colors)
    "CO2": "#E69F00",     # Orange"
    "H2":  "#332288",   # Indigo
    "H2O": "#56B4E9",   # Sky blue
    "N2":  "#882255",   # Purple 
    
    # Other major species
    "H": "#D55E00",   # Vermillion (orange-red)
    "O": "#009E73",     # Bluish green
    "C": "#F0E442",     # Yellow
    "CH4":"#117733",   # Green
    "CO": "#AA4499",    # Mauve
    "O2": "#44AA99",    # Teal
    "N": "#999933",     # Olive
    "NH3": "#CC79A7",    # Reddish purple
    "H2S": "#88CCEE",   # Cyan
    "S": "#DDCC77",     # Sand
    "S2": "#661100",    # Dark red
    "SO2": "#332288",    # Blue
}

mass_labels = {
    "1_M_earth": r"$M_p = 1\,M_{\oplus}$",
    "10_M_earth": r"$M_p = 10\,M_{\oplus}$"
}

ratios = []

for case, group in df.groupby(case_cols):
    group = group[group["Mdot_kg_s"] > 0].sort_values("Mdot_kg_s", ascending=False)

    if len(group) < 2:
        continue

    dominant = group.iloc[0]
    second = group.iloc[1]

    ratio = dominant["Mdot_kg_s"] / second["Mdot_kg_s"]

    ratios.append({
        "file": case[0],
        "T_inf": case[1],
        "atmosphere_type": case[2],
        "mass_case": case[3],
        "flux_case": case[4],
        "dominant_species": dominant["species"],
        "second_species": second["species"],
        "dominant_Mdot": dominant["Mdot_kg_s"],
        "second_Mdot": second["Mdot_kg_s"],
        "dominant_to_second_ratio": ratio,
        "log10_ratio": np.log10(ratio),
    })

ratio_df = pd.DataFrame(ratios)
ratio_df.to_csv("dominant_to_second_species_ratios.csv", index=False)

plt.figure(figsize=(7, 5))
plt.hist(ratio_df["log10_ratio"], bins=20, edgecolor="black")

plt.xlabel(r"$\log_{10}(\dot{M}_{\rm dominant}/\dot{M}_{\rm second})$", fontsize=13)
plt.ylabel("Number of cases", fontsize=13)
plt.title("Dominance of the largest escaping species", fontsize=14)
plt.tight_layout()
plt.savefig("Plots/PROTEUS plots/hist_dominant_to_second_ratio.png", dpi=300, bbox_inches="tight")
plt.show()



full = pd.read_csv("proteus_jeans_escape_results.csv")

ratio_rows = []

group_cols = [
    "atmosphere_type",
    "mass_case",
    "flux_case",
    "T_inf",
]

for keys, group in full.groupby(group_cols):

    group = group.sort_values("Mdot_kg_s", ascending=False)

    if len(group) < 2:
        continue

    dominant = group.iloc[0]
    second = group.iloc[1]

    ratio = dominant["Mdot_kg_s"] / second["Mdot_kg_s"]

    ratio_rows.append({
        "atmosphere_type": keys[0],
        "mass_case": keys[1],
        "flux_case": keys[2],
        "T_inf": keys[3],

        "dominant_species": dominant["species"],
        "second_species": second["species"],

        "ratio": ratio,
        "log_ratio": np.log10(ratio),
    })

ratio_df = pd.DataFrame(ratio_rows)

from matplotlib.lines import Line2D

flux_legend = [
    Line2D([0], [0], marker="o", color="black", linestyle="None",
           label=r"1 $F_{\oplus}$"),
    Line2D([0], [0], marker="^", color="black", linestyle="None",
           label=r"1000 $F_{\oplus}$"),
]

atm_handles = [
    Line2D([0], [0], color=species_colors[atm], lw=2, label=atm)
    for atm in atm_archetype
]


fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for idx, mas in enumerate(mass):

    ax = axes[idx]

    for atm in atm_archetype:

        for inst in instellation:

            s = ratio_df[
                (ratio_df["atmosphere_type"] == atm) &
                (ratio_df["mass_case"] == mas) &
                (ratio_df["flux_case"] == inst)
            ].sort_values("T_inf")

            if s.empty:
                continue

            ax.scatter(
                s["T_inf"],
                s["log_ratio"],
                color=species_colors[atm],
                marker="o" if inst == "1_F_earth" else "^",
                s=70,
            )

            ax.scatter(
                s["T_inf"],
                s["log_ratio"],
                color=species_colors[atm],
                alpha=0.5,
            )

            for _, row in s.iterrows():
                label = f"{row['dominant_species']}/{row['second_species']}"

                ax.annotate(
                    label,
                    (row["T_inf"], row["log_ratio"]),
                    fontsize=7,
                    xytext=(0, 4),
                    textcoords="offset points",
                    ha="center",
                )

    ax.set_xlabel(r"$T_{\infty}$ [K]")
    ax.set_title(mass_labels[mas])
    first_legend = axes[1].legend(handles=atm_handles, title="Atmosphere", bbox_to_anchor=(0.80, 1))
    axes[1].add_artist(first_legend)
    axes[1].legend(handles=flux_legend, title="Instellation")
    first_legend = axes[0].legend(handles=atm_handles, title="Atmosphere", bbox_to_anchor=(0.80, 1))
    axes[0].add_artist(first_legend)
    axes[0].legend(handles=flux_legend, title="Instellation")

axes[0].set_ylabel(
    r"$\log_{10}(\dot M_{\rm dominant}/\dot M_{\rm second})$"
)

plt.suptitle("Dominance of largest escaping species")
plt.tight_layout()
plt.savefig(
    "Plots/PROTEUS plots/dominant_second_ratio_vs_Tinf.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()