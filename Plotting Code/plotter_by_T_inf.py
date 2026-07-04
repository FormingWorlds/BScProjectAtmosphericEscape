import os
import matplotlib.pyplot as plt
import scienceplots
#plt.style.use('science') 
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

outdir = "Plots/PROTEUS_Bates_sensitivity/Dissociation" #change either Dissociation or Constant VMR 
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

atm_handles = [
    Line2D([0], [0], color=species_colors[atm], lw=2, label=atm)
    for atm in atm_archetype
]

invalid_handle = Line2D(
    [0], [0],
    marker="x",
    color="black",
    linestyle="None",
    markersize=8,
    markeredgewidth=2,
    label=r"Outside Jeans Validity ($\lambda_J < 1.5$)"
)


all_handles = (
    atm_handles
    + flux_handles
    + [invalid_handle]
)

df = pd.read_csv("Outputs/proteus_jeans_case_summary_dissociation.csv") #change if dissociation

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

                valid = s[
                    (s["jeans_valid"])
                    &
                    (~s["unphysical_extension"])
                ]

                hydro = s[
                    (~s["jeans_valid"])
                    &
                    (~s["unphysical_extension"])
                ]


                if s.empty:
                    continue

                #ax.plot(
                    #valid["T_inf"],
                    #valid[ycol],
                    #color=species_colors[atm],
                    #linestyle="-",
                    #linewidth=1.5,
                    #alpha=0.7,
                #)

                ax.scatter(
                    valid["T_inf"],
                    valid[ycol],
                    color=species_colors[atm],
                    marker=flux_markers[flux],
                    s=65,
                )

                ax.scatter(
                    hydro["T_inf"],
                    hydro[ycol],
                    color=species_colors[atm],
                    marker="x",
                    label="Hydrodynamic onset",
                    s=90,
                    linewidths=2,
                )

        if ylog:
            ax.set_yscale("log")

        ax.set_xlabel(r"$T_{\infty}$ [K]", fontsize=13)
        ax.set_title(mass_labels[mas], fontsize=14)
        ax.tick_params(axis="both", which="major", labelsize=12)

    axes[0].set_ylabel(ylabel, fontsize=13)

    fig.legend(
        handles=all_handles,
        loc="lower center",
        title="Legend",
        bbox_to_anchor=(0.5, -0.02),
        ncol=len(all_handles),
        frameon=False,
        fontsize=12,
    )
    #axes[1].legend(handles=flux_handles + [invalid_handle], title="Markers", loc="lower right")

    #axes[0].legend(handles=flux_handles, title="Instellation", loc="lower right")

    plt.tight_layout(rect=[0, 0.08, 1, 1])
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


#5. Weighted mass loss 

df["log10_weighted_mass_loss"] = np.log10(
    df["weighted_mass_loss_kg_s"].replace(0, np.nan)
)

plot_vs_Tinf(
    ycol="log10_weighted_mass_loss",
    ylabel=r"$\log_{10}(\dot{M}_{\rm weighted})$ [kg/s]",
    filename="weighted_mass_loss_vs_Tinf.png",
    ylog=False,
)


# 6. Dominant-species lambda_J vs T_inf

full = pd.read_csv("Outputs/proteus_jeans_escape_results_dissociation.csv")

# Make sure column exists for older files
if "unphysical_extension" not in full.columns:
    full["unphysical_extension"] = False

# Remove placeholder rows before choosing dominant species
physical_full = full[~full["unphysical_extension"]].copy()

dominant = (
    physical_full.sort_values("Mdot_kg_s", ascending=False)
    .groupby(["atmosphere_type", "mass_case", "flux_case", "T_inf"], as_index=False)
    .first()
)

# Add unphysical rows separately from full file
unphysical_cases = full[full["unphysical_extension"]].copy()

fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for idx, mas in enumerate(mass):
    ax = axes[idx]

    for atm in atm_archetype:
        for inst in instellation:

            s = dominant[
                (dominant["atmosphere_type"] == atm) &
                (dominant["mass_case"] == mas) &
                (dominant["flux_case"] == inst)
            ].sort_values("T_inf")


            if s.empty:
                continue

            valid = s[s["jeans_valid"]]
            hydro = s[~s["jeans_valid"]]

            ax.scatter(
                valid["T_inf"],
                valid["lambda_j"],
                color=species_colors.get(atm),
                marker=flux_markers[inst],
                s=80,
            )

            ax.scatter(
                hydro["T_inf"],
                hydro["lambda_j"],
                color=species_colors.get(atm),
                marker="x",
                s=90,
                linewidths=2,
            )

            for _, row in s.iterrows():
                ax.annotate(
                    row["species"],
                    (row["T_inf"], row["lambda_j"]),
                    fontsize=8,
                    xytext=(0, 3.5),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                )

    ax.axhline(1.5, color="grey", linestyle="--", linewidth=1)
    ax.axhline(10, color="grey", linestyle=":", linewidth=1)

    ax.set_yscale("log")
    ax.set_xlabel(r"$T_{\infty}$ [K]", fontsize=13)
    ax.set_title(mass_labels[mas], fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)

axes[0].set_ylabel(r"Dominant species Jeans parameter $\lambda_J$", fontsize=13)

fig.legend(
    handles=all_handles,
    loc="lower center",
    title="Legend",
    bbox_to_anchor=(0.5, -0.02),
    ncol=len(all_handles),
    frameon=False,
    fontsize=12,
)

plt.suptitle(r"Dominant-species Jeans parameter vs $T_{\infty}$", fontsize=15)
plt.tight_layout(rect=[0, 0.08, 1, 1])
plt.savefig("Plots/PROTEUS_Bates_sensitivity/Dissociation/dominant_lambda_vs_Tinf.png", dpi=300, bbox_inches="tight")
plt.close()


# 7. Mass loss vs dominant-species lambda_J

fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for idx, mas in enumerate(mass):
    ax = axes[idx]

    for atm in atm_archetype:
        for inst in instellation:

            s = dominant[
                (dominant["atmosphere_type"] == atm) &
                (dominant["mass_case"] == mas) &
                (dominant["flux_case"] == inst)
            ]

            if s.empty:
                continue

            valid = s[s["jeans_valid"]]
            hydro = s[~s["jeans_valid"]]

            ax.scatter(
                valid["lambda_j"],
                np.log10(valid["Mdot_kg_s"].replace(0, np.nan)),
                color=species_colors.get(atm),
                marker=flux_markers[inst],
                s=80,
            )

            ax.scatter(
                hydro["lambda_j"],
                np.log10(hydro["Mdot_kg_s"].replace(0, np.nan)),
                color=species_colors.get(atm),
                marker="x",
                s=90,
                linewidths=2,
            )

    ax.axvline(1.5, color="grey", linestyle="--", linewidth=1)
    ax.axvline(10, color="grey", linestyle=":", linewidth=1)

    ax.set_xscale("log")
    ax.invert_xaxis()
    ax.set_xlabel(r"Dominant species $\lambda_J$", fontsize=13)
    ax.set_title(mass_labels[mas], fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)

axes[0].set_ylabel(r"$\log_{10}(\dot{M}_{\rm dominant})$ [kg/s]", fontsize=13)

fig.legend(
    handles=all_handles,
    loc="lower center",
    title="Legend",
    bbox_to_anchor=(0.5, -0.02),
    ncol=len(all_handles),
    frameon=False,
    fontsize=12,
)

plt.suptitle(r"Dominant species mass loss versus Jeans parameter", fontsize=15)
plt.tight_layout(rect=[0, 0.08, 1, 1])
plt.savefig("Plots/PROTEUS_Bates_sensitivity/Dissociation/dominant_mass_loss_vs_lambda.png", dpi=300, bbox_inches="tight")
plt.close()

# 8. Hydrogen lambda_J vs T_inf

hydrogen = physical_full[physical_full["species"] == "H"].copy()
unphysical_h = full[full["unphysical_extension"]].copy()

fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for idx, mas in enumerate(mass):
    ax = axes[idx]

    for atm in atm_archetype:
        for inst in instellation:

            s = hydrogen[
                (hydrogen["atmosphere_type"] == atm) &
                (hydrogen["mass_case"] == mas) &
                (hydrogen["flux_case"] == inst)
            ].sort_values("T_inf")

            if s.empty:
                continue

            valid = s[s["jeans_valid"]]
            hydro = s[~s["jeans_valid"]]

            ax.plot(
                valid["T_inf"],
                valid["lambda_j"],
                color=species_colors.get(atm),
                linewidth=1.8,
                alpha=0.8,
            )

            ax.scatter(
                valid["T_inf"],
                valid["lambda_j"],
                color=species_colors.get(atm),
                marker=flux_markers[inst],
                s=60,
            )

            ax.scatter(
                hydro["T_inf"],
                hydro["lambda_j"],
                color=species_colors.get(atm),
                marker="x",
                s=90,
                linewidths=2,
            )


    ax.axhline(1.5, color="grey", linestyle="--", linewidth=1)
    ax.axhline(10, color="grey", linestyle=":", linewidth=1)

    ax.set_yscale("log")
    ax.set_xlabel(r"$T_{\infty}$ [K]", fontsize=13)
    ax.set_title(mass_labels[mas], fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)

axes[0].set_ylabel(r"Hydrogen Jeans parameter $\lambda_{J,\mathrm{H}}$", fontsize=13)

fig.legend(
    handles=all_handles,
    loc="lower center",
    title="Legend",
    bbox_to_anchor=(0.5, -0.02),
    ncol=len(all_handles),
    frameon=False,
    fontsize=12,
)

plt.suptitle(r"Hydrogen Jeans parameter vs $T_{\infty}$", fontsize=15)
plt.tight_layout(rect=[0, 0.08, 1, 1])
plt.savefig("Plots/PROTEUS_Bates_sensitivity/Dissociation/hydrogen_lambda_vs_Tinf.png", dpi=300, bbox_inches="tight")
plt.close()