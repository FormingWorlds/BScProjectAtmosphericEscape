import os
import matplotlib.pyplot as plt
import scienceplots
#plt.style.use('science') 
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

outdir = "Plots/PROTEUS_Bates_sensitivity"
os.makedirs(outdir, exist_ok=True)

atm_archetype = ["CO2", "N2", "H2O", "H2"]
instellation = ["1_F_earth", "1000_F_earth"]
mass = ["1_M_earth", "10_M_earth"]

mass_labels = {
    "1_M_earth": r"$M_p = 1\,M_{\oplus}$",
    "10_M_earth": r"$M_p = 10\,M_{\oplus}$"
}

flux_markers = {
    "1_F_earth": "o",
    "1000_F_earth": "^",
}

flux_labels = {
    "1_F_earth": r"1 $F_{\oplus}$",
    "1000_F_earth": r"1000 $F_{\oplus}$",
}

flux_handles = [
    Line2D([0], [0], marker=flux_markers[f], color="black",
           linestyle="None", label=flux_labels[f])
    for f in instellation
]



# Define consistent colorblind-friendly colors for each species
# Using Paul Tol's colorblind-safe palette
species_colors = {
    # Atmosphere archetypes (distinctive colors)
    "CO2": "#D55E00",   # Vermillion (orange-red)
    "H2":  "#332288",   # Indigo
    "H2O": "#56B4E9",   # Sky blue
    "N2":  "#117733",   # Green
    
    # Other major species
    "H": "#E69F00",     # Orange
    "O": "#009E73",     # Bluish green
    "C": "#F0E442",     # Yellow
    "CH4": "#882255",   # Purple
    "CO": "#AA4499",    # Mauve
    "O2": "#44AA99",    # Teal
    "N": "#999933",     # Olive
    "NH3": "#CC79A7",    # Reddish purple
    "H2S": "#88CCEE",   # Cyan
    "S": "#DDCC77",     # Sand
    "S2": "#661100",    # Dark red
    "SO2": "#332288",    # Blue
}

atm_handles = [
    Line2D([0], [0], color=species_colors[atm], lw=2, label=atm)
    for atm in atm_archetype
]


df = pd.read_csv("proteus_jeans_case_summary.csv")
df["log10_weighted_mass_loss"] = np.log10(df["weighted_mass_loss_kg_s"])

def plot_vs_Tinf(ycol, ylabel, filename, ylog=False):

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    for idx, mas in enumerate(mass):
        ax = axes[idx]

        for atm in atm_archetype:
            for flux in instellation:
                s = df[
                    (df["mass_case"] == mas) &
                    (df["atmosphere_type"] == atm) &
                    (df["flux_case"] == flux)
                ].sort_values("T_inf")

                if s.empty:
                    continue

                ax.scatter(
                    s["T_inf"],
                    s[ycol],
                    color=species_colors[atm],
                    marker=flux_markers[flux],
                    linestyle="-",
                    linewidth=1.8,
                )

        if ylog:
            ax.set_yscale("log")

        ax.set_xlabel(r"$T_{\infty}$ [K]", fontsize=13)
        ax.set_title(mass_labels[mas], fontsize=14)
        ax.tick_params(axis="both", which="major", labelsize=12)

    axes[0].set_ylabel(ylabel, fontsize=13)

    fig.legend(handles=atm_handles, title="Atmosphere", bbox_to_anchor=(1.10, 1), loc="outside upper right")

    axes[1].legend(handles=flux_handles, title="Instellation", loc="lower right")

    axes[0].legend(handles=flux_handles, title="Instellation", loc="lower right")

    plt.tight_layout()
    plt.savefig(os.path.join(outdir, filename), dpi=300, bbox_inches="tight")
    plt.close()


# 1. Main sensitivity plot
plot_vs_Tinf(
    ycol="log10_weighted_mass_loss",
    ylabel=r"$\log_{10}(\dot{M}_{\rm weighted})$ [kg/s]",
    filename="weighted_mass_loss_vs_Tinf.png",
    ylog=False,
)

# 2. Exobase altitude sensitivity
plot_vs_Tinf(
    ycol="exobase_altitude_km",
    ylabel="Exobase altitude [km]",
    filename="exobase_altitude_vs_Tinf.png",
    ylog=True,
)

# 3. Exobase temperature sanity check
plot_vs_Tinf(
    ycol="T_exo_K",
    ylabel=r"$T_{\rm exo}$ [K]",
    filename="Texo_vs_Tinf.png",
    ylog=False,
)

# 4. Dominant species mass loss sensitivity
df["log10_dominant_mass_loss"] = np.log10(
    df["dominant_species_Mdot_kg_s"].replace(0, np.nan)
)

plot_vs_Tinf(
    ycol="log10_dominant_mass_loss",
    ylabel=r"$\log_{10}(\dot{M}_{\rm dominant})$ [kg/s]",
    filename="dominant_mass_loss_vs_Tinf.png",
    ylog=False,
)