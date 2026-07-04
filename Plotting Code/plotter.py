import matplotlib.pyplot as plt
import scienceplots
#plt.style.use('science') 
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

atm_archetype = ["CO2", "N2", "H2O", "H2"]
instellation = ["1_F_earth", "1000_F_earth"]
mass = ["1_M_earth", "10_M_earth"]

mass_labels = {
    "1_M_earth": r"$M_p = 1\,M_{\oplus}$",
    "10_M_earth": r"$M_p = 10\,M_{\oplus}$"
}

instellation_labels = {
    "1_F_earth": r"1 $F_{\oplus}$",
    "1000_F_earth": r"1000 $F_{\oplus}$"
}

flux_legend = [
    Line2D([0], [0], marker="o", color="black", linestyle="None", label=r"1 $F_{\oplus}$"),
    Line2D([0], [0], marker="^", color="black", linestyle="None", label=r"1000 $F_{\oplus}$"),
]


# Define consistent colorblind-friendly colors for each species
# Using Paul Tol's colorblind-safe palette
species_colors = {
    # Atmosphere archetypes (distinctive colors)
    "CO2": "#D55E00",   # Vermillion (orange-red)
    "H2": "#0072B2",    # Blue
    "H2O": "#56B4E9",   # Sky blue
    "N2": "#CC79A7",    # Reddish purple
    
    # Other major species
    "H": "#E69F00",     # Orange
    "O": "#009E73",     # Bluish green
    "C": "#F0E442",     # Yellow
    "CH4": "#882255",   # Purple
    "CO": "#AA4499",    # Mauve
    "O2": "#44AA99",    # Teal
    "N": "#999933",     # Olive
    "NH3": "#117733",   # Green
    "H2S": "#88CCEE",   # Cyan
    "S": "#DDCC77",     # Sand
    "S2": "#661100",    # Dark red
    "SO2": "#332288"    # Indigo
}

atm_handles = [
    Line2D([0], [0], marker="o", color=species_colors[atm], linestyle="None", label=atm)
    for atm in atm_archetype
]

df = pd.read_csv("Outputs/proteus_jeans_case_summary.csv")
df["log10_weighted_mass_loss"] = np.log10(df["weighted_mass_loss_kg_s"])

#print(df[["atmosphere_type", "mass_case", "flux_case"]])

# 1. Weighted mass loss per atmosphere

fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for idx, mas in enumerate(mass):
    ax = axes[idx]

    for atm in atm_archetype:
        for inst in instellation:
            subset = df[
                (df["atmosphere_type"] == atm) &
                (df["mass_case"] == mas) &
                (df["flux_case"] == inst)
            ]

            if subset.empty:
                continue

            y = (subset["log10_weighted_mass_loss"].iloc[0])

            ax.scatter(
                atm,
                y,
                color=species_colors.get(atm),
                marker="o" if inst == "1_F_earth" else "^",
                s=80,
            )

    
    ax.set_xlabel("Atmosphere type", fontsize=13)
    ax.set_ylabel(r"$\log_{10}(\dot{M}_{weighted})$ [kg/s]")
    ax.set_title(mass_labels[mas], fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)
    first_legend = axes[1].legend(handles=atm_handles, title="Atmosphere", loc="upper left")
    axes[1].add_artist(first_legend)
    first_legend = axes[0].legend(handles=atm_handles, title="Atmosphere", loc="upper left")
    axes[0].add_artist(first_legend)

    axes[1].legend(handles=flux_legend, title="Instellation", fontsize=11, loc="lower right")
    axes[0].legend(handles=flux_legend, title="Instellation", fontsize=11, loc ="lower right")

plt.suptitle("Weighted Jeans escape by atmosphere type", fontsize=15)
plt.tight_layout()
plt.savefig("Plots/PROTEUS plots/weighted_mass_loss_by_atmosphere.png", dpi=300, bbox_inches="tight")
plt.close()


# 2. Dominant species mass loss per atmosphere

fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for idx, mas in enumerate(mass):
    ax = axes[idx]

    for atm in atm_archetype:
        for inst in instellation:
            subset = df[
                (df["atmosphere_type"] == atm) &
                (df["mass_case"] == mas) &
                (df["flux_case"] == inst)
            ]

            if subset.empty:
                continue

            y = np.log10(subset["dominant_species_Mdot_kg_s"].iloc[0])
            dom_sp = subset["dominant_escaping_species"].iloc[0]

            ax.scatter(
                atm,
                y,
                color=species_colors.get(dom_sp, "black"),
                marker="o" if inst == "1_F_earth" else "^",
                s=80,
            )

            ax.annotate(
                dom_sp,
                (atm,y),
                fontsize=9,
                xytext=(0, 3.5),               # Offsets the text by 0 points horizontally and 5 points vertically
                textcoords='offset points',
                ha="center",
                va="bottom",
            )

    ax.set_xlabel("Atmosphere type", fontsize=13)
    ax.set_ylabel(r"Dominant species $\log_{10}(\dot{M}_{dot})$ [kg/s]")
    ax.set_title(mass_labels[mas], fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)

    axes[1].legend(handles=flux_legend, title="Instellation", fontsize=11)
    axes[0].legend(handles=flux_legend, title="Instellation", fontsize=11)

plt.suptitle("Dominant escaping species", fontsize=15)
plt.tight_layout()
plt.savefig("Plots/PROTEUS plots/dominant_species_mass_loss_by_atmosphere.png", dpi=300, bbox_inches="tight")
plt.close()



# 3. Exobase altitude per atmosphere

fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for idx, mas in enumerate(mass):
    ax = axes[idx]

    for atm in atm_archetype:
        for inst in instellation:
            subset = df[
                (df["atmosphere_type"] == atm) &
                (df["mass_case"] == mas) &
                (df["flux_case"] == inst)
            ]

            if subset.empty:
                continue

            y = subset["exobase_altitude_km"].iloc[0]

            ax.scatter(
                atm,
                y,
                color=species_colors.get(atm),
                marker="o" if inst == "1_F_earth" else "^",
                s=80,
            )
    ax.set_yscale('log')
    ax.set_xlabel("Atmosphere type", fontsize=13)
    ax.set_title(mass_labels[mas], fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)

    axes[0].set_ylabel("Exobase altitude [km]", fontsize=13)
    axes[1].set_ylabel("Exobase altitude [km]", fontsize=13)
    first_legend = axes[1].legend(handles=atm_handles, title="Atmosphere", loc="upper left")
    axes[1].add_artist(first_legend)
    first_legend = axes[0].legend(handles=atm_handles, title="Atmosphere", loc="upper left")
    axes[0].add_artist(first_legend)

    axes[1].legend(handles=flux_legend, title="Instellation", loc="lower right", fontsize=11)
    axes[0].legend(handles=flux_legend, title="Instellation", loc="upper right", fontsize=11)

plt.suptitle("Exobase altitude by atmosphere type", fontsize=15)
plt.tight_layout()
plt.savefig("Plots/PROTEUS plots/exobase_altitude_by_atmosphere.png", dpi=300, bbox_inches="tight")
plt.close()


# 4. Temperature vs weighted mass loss

fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for idx, mas in enumerate(mass):
    ax = axes[idx]

    for atm in atm_archetype:
        subset = df[df["atmosphere_type"].eq(atm) & df["mass_case"].eq(mas)]

        if subset.empty:
            continue

        for inst in instellation:
            s = subset[subset["flux_case"] == inst]

            ax.scatter(
                s["T_exo_K"],
                np.log10(s["weighted_mass_loss_kg_s"]),
                color=species_colors.get(atm),
                marker="o" if inst == "1_F_earth" else "^",
                s=80,
                label=atm if inst == "1_F_earth" else None,
            )

    ax.set_xlabel("Exobase temperature [K]", fontsize=13)
    ax.set_title(mass_labels[mas], fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)

    axes[0].set_ylabel(r"$\log_{10}(\dot{M}_{weighted})$ [kg/s]", fontsize=13)
    axes[1].set_ylabel(r"$\log_{10}(\dot{M}_{weighted})$ [kg/s]", fontsize=13)

    first_legend = axes[1].legend(handles=atm_handles, title="Atmosphere", bbox_to_anchor=(1.0, 0.55))
    axes[1].add_artist(first_legend)
    axes[1].legend(handles=flux_legend, title="Instellation", loc="lower right")
    first_legend = axes[0].legend(handles=atm_handles, title="Atmosphere",  bbox_to_anchor=(1.0, 0.55))
    axes[0].add_artist(first_legend)
    axes[0].legend(handles=flux_legend, title="Instellation", loc="lower right")

plt.suptitle("Mass loss sensitivity to exobase temperature", fontsize=15)
plt.tight_layout()
plt.savefig("Plots/PROTEUS plots/mass_loss_vs_exobase_temperature.png", dpi=300, bbox_inches="tight")
plt.close()



# 5. Exobase altitude vs weighted mass loss

fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for idx, mas in enumerate(mass):
    ax = axes[idx]

    for atm in atm_archetype:
        subset = df[df["atmosphere_type"].eq(atm) & df["mass_case"].eq(mas)]

        if subset.empty:
            continue

        for inst in instellation:
            s = subset[subset["flux_case"] == inst]

            ax.scatter(
                s["exobase_altitude_km"],
                np.log10(s["weighted_mass_loss_kg_s"]),
                color=species_colors.get(atm),
                marker="o" if inst == "1_F_earth" else "^",
                s=80,
            )

    ax.set_xscale("log")
    ax.set_xlabel("Exobase altitude [km]", fontsize=13)
    ax.set_title(mass_labels[mas], fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)

    axes[0].set_ylabel(r"$\log_{10}(\dot{M}_{weighted})$ [kg/s]", fontsize=13)
    axes[1].set_ylabel(r"$\log_{10}(\dot{M}_{weighted})$ [kg/s]", fontsize=13)

    first_legend = axes[1].legend(handles=atm_handles, title="Atmosphere", bbox_to_anchor=(1.0, 0.55))
    axes[1].add_artist(first_legend)
    axes[1].legend(handles=flux_legend, title="Instellation", loc="lower right")
    first_legend = axes[0].legend(handles=atm_handles, title="Atmosphere",bbox_to_anchor=(1.0, 0.55))
    axes[0].add_artist(first_legend)
    axes[0].legend(handles=flux_legend, title="Instellation", loc="lower right")

plt.suptitle("Mass loss sensitivity to exobase altitude", fontsize=15)
plt.tight_layout()
plt.savefig("Plots/PROTEUS plots/mass_loss_vs_exobase_altitude.png", dpi=300, bbox_inches="tight")
plt.close()
