import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.constants import G, pi

outdir = "Plots/Predictions"
os.makedirs(outdir, exist_ok=True)

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

# Load lower and upper bound cases
constant = pd.read_csv("proteus_jeans_case_summary.csv")
diss = pd.read_csv("proteus_jeans_case_summary_dissociation.csv")

constant["case"] = "Constant VMR"
diss["case"] = "Full dissociation"

case_styles = {
    "Constant VMR": "-",
    "Full dissociation": "--",
}

df = pd.concat([constant, diss], ignore_index=True)

df["log10_weighted_mass_loss"] = np.log10(
    df["weighted_mass_loss_kg_s"].replace(0, np.nan)
)

invalid_handle = Line2D(
    [0], [0],
    marker="x",
    color="black",
    linestyle="None",
    markersize=8,
    markeredgewidth=2,
    label=r"Hydrodynamic onset ($\lambda_J < 1.5$)"
)

case_handles = [
    Line2D([0], [0], color="black", linestyle="-", label="Constant VMR"),
    Line2D([0], [0], color="black", linestyle="--", label="Full dissociation"),
]

flux_handles = [
    Line2D([0], [0], marker=flux_markers[f], color="black",
           linestyle="None", label=flux_labels[f])
    for f in flux_cases
]

atm_handles = [
    Line2D([0], [0], color=species_colors[a], lw=2, label=a)
    for a in atm_archetype
]



fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for idx, mass in enumerate(mass_cases):
    ax = axes[idx]

    for atm in atm_archetype:
        for flux in flux_cases:
            for case, linestyle in [
                ("Constant VMR", "-"),
                ("Full dissociation", "--"),
            ]:
                s = df[
                    (df["mass_case"] == mass) &
                    (df["atmosphere_type"] == atm) &
                    (df["flux_case"] == flux) &
                    (df["case"] == case)
                ].sort_values("T_inf")

                if s.empty:
                    continue

                valid = s[s["jeans_valid"]]
                invalid = s[~s["jeans_valid"]]

                ax.plot(
                    valid["T_inf"],
                    valid["log10_weighted_mass_loss"],
                    color=species_colors[atm],
                    linestyle=linestyle,
                    linewidth=1.7,
                    alpha=0.85,
                )

                ax.scatter(
                    valid["T_inf"],
                    valid["log10_weighted_mass_loss"],
                    color=species_colors[atm],
                    marker=flux_markers[flux],
                    s=50,
                )

                ax.scatter(
                    invalid["T_inf"],
                    invalid["log10_weighted_mass_loss"],
                    color=species_colors[atm],
                    marker="x",
                    s=80,
                    linewidths=2,
                )

    ax.set_xlabel(r"$T_{\infty}$ [K]", fontsize=13)
    ax.set_title(mass_labels[mass], fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)

axes[0].set_ylabel(r"$\log_{10}(\dot{M}_{\rm weighted})$ [kg/s]", fontsize=13)

all_handles = atm_handles + case_handles + flux_handles + [invalid_handle]

fig.legend(
    handles=all_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.12),
    ncol=4,
    frameon=False,
    fontsize=9,
)

#plt.suptitle("Predicted Jeans escape bracket", fontsize=15)
plt.tight_layout(rect=[0, 0.12, 1, 1])
plt.savefig(os.path.join(outdir, "escape_bracket_vs_Tinf.png"), dpi=300, bbox_inches="tight")
plt.close()


fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for idx, mass in enumerate(mass_cases):
    ax = axes[idx]

    for atm in atm_archetype:
        for flux in flux_cases:
            for case, linestyle in case_styles.items():

                s = df[
                    (df["mass_case"] == mass) &
                    (df["atmosphere_type"] == atm) &
                    (df["flux_case"] == flux) &
                    (df["case"] == case)
                ].sort_values("T_inf")

                if s.empty:
                    continue

                valid = s[s["jeans_valid"]]
                invalid = s[~s["jeans_valid"]]

                ax.plot(
                    s["T_inf"],
                    s["dominant_lambda_j"],
                    color=species_colors[atm],
                    linestyle=linestyle,
                    linewidth=1.6,
                    alpha=0.85,
                )

                ax.scatter(
                    valid["T_inf"],
                    valid["dominant_lambda_j"],
                    color=species_colors[atm],
                    marker=flux_markers[flux],
                    s=45,
                )

                ax.scatter(
                    invalid["T_inf"],
                    invalid["dominant_lambda_j"],
                    color=species_colors[atm],
                    marker="x",
                    s=80,
                    linewidths=2,
                )

    ax.axhline(1.5, color="grey", linestyle="--", linewidth=1)
    ax.axhline(10, color="grey", linestyle=":", linewidth=1)

    ax.set_yscale("log")
    ax.set_xlabel(r"$T_{\infty}$ [K]", fontsize=13)
    ax.set_title(mass_labels[mass], fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)

axes[0].set_ylabel(r"Dominant species Jeans parameter $\lambda_J$", fontsize=13)

fig.legend(
    handles=all_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.2),
    ncol=3,
    frameon=False,
    fontsize=12,
)

plt.suptitle("Transition from Jeans to non-Jeans escape", fontsize=15)
plt.tight_layout(rect=[0, 0.16, 1, 1])
plt.savefig(os.path.join(outdir, "dominant_lambda_transition_vs_Tinf.png"), dpi=300, bbox_inches="tight")
plt.close()


# Load planet radii/masses from bulk files
def read_bulk(comp, mass_case, flux_case):
    path = (
        f"PROTEUS data/{comp}_atmospheres/"
        f"planet_bulk_properties_{comp}_atmospheres_{mass_case}.csv"
    )

    bulk = pd.read_csv(path, sep="\t")
    row = bulk[bulk["Case"] == flux_case].iloc[0]

    return row["R_int [m]"], row["M_planet [kg]"]


P_SURF = 1e5          # Pa, 1-bar equivalent atmosphere
SECONDS_PER_YEAR = 365.25 * 24 * 3600

def one_bar_atmosphere_mass(R_p, M_p):
    g = G * M_p / R_p**2
    return 4 * np.pi * R_p**2 * P_SURF / g


# Add 1-bar equivalent atmosphere mass and lifetime
M_atm_values = []

for _, row in df.iterrows():
    R_p, M_p = read_bulk(
        row["atmosphere_type"],
        row["mass_case"],
        row["flux_case"],
    )

    M_atm = one_bar_atmosphere_mass(R_p, M_p)
    M_atm_values.append(M_atm)

df["M_atm_1bar_kg"] = M_atm_values
df["lifetime_s"] = df["M_atm_1bar_kg"] / df["weighted_mass_loss_kg_s"]
df["lifetime_yr"] = df["lifetime_s"] / SECONDS_PER_YEAR
df["log10_lifetime_yr"] = np.log10(df["lifetime_yr"].replace([np.inf, 0], np.nan))

df.to_csv("prediction_lifetimes_1bar.csv", index=False)

# Legend handles
atm_handles = [
    Line2D([0], [0], color=species_colors[a], lw=2, label=a)
    for a in atm_archetype
]

reference_handles = [
    Line2D([0],[0], linestyle=":", color="gray", label="1 Myr"),
    Line2D([0],[0], linestyle="--", color="gray", label="1 Gyr"),
    Line2D([0],[0], linestyle="-.", color="gray", label="Earth age"),
]

all_handles = atm_handles + flux_handles + [invalid_handle] + reference_handles

fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for idx, mass in enumerate(mass_cases):
    ax = axes[idx]

    for atm in atm_archetype:
        for flux in flux_cases:
            for case, linestyle in case_styles.items():

                s = df[
                    (df["mass_case"] == mass) &
                    (df["atmosphere_type"] == atm) &
                    (df["flux_case"] == flux) &
                    (df["case"] == case)
                ].sort_values("T_inf")

                if s.empty:
                    continue

                valid = s[s["jeans_valid"]]
                invalid = s[~s["jeans_valid"]]


                ax.scatter(
                    valid["T_inf"],
                    valid["log10_lifetime_yr"],
                    color=species_colors[atm],
                    marker=flux_markers[flux],
                    s=45,
                )

                ax.scatter(
                    invalid["T_inf"],
                    invalid["log10_lifetime_yr"],
                    color=species_colors[atm],
                    marker="x",
                    s=80,
                    linewidths=2,
                )

    # Useful reference timescales
    ax.axhline(6, color="grey", linestyle=":", linewidth=1, label="1Myr")       # 1 Myr
    ax.axhline(9, color="grey", linestyle="--", linewidth=1, label="1Gyr")      # 1 Gyr
    ax.axhline(np.log10(4.5e9), color="grey", linestyle="-.", linewidth=1, label="Earth age")  # Earth age

    ax.set_xlabel(r"$T_{\infty}$ [K]", fontsize=13)
    ax.set_title(mass_labels[mass], fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)

axes[0].set_ylabel(
    r"$\log_{10}(\tau_{\rm loss,\,1bar})$ [yr]",
    fontsize=13,
)

fig.legend(
    handles=all_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.03),
    ncol=3,
    frameon=False,
    fontsize=12,
)

plt.suptitle("Predicted lifetime of a 1-bar equivalent atmosphere", fontsize=15)
plt.tight_layout(rect=[0, 0.16, 1, 1])
plt.savefig(os.path.join(outdir, "lifetime_1bar_vs_Tinf.png"), dpi=300, bbox_inches="tight")
plt.close()



const = pd.read_csv("proteus_jeans_case_summary.csv")
diss = pd.read_csv("proteus_jeans_case_summary_dissociation.csv")

merge_cols = [
    "file",
    "T_inf",
    "atmosphere_type",
    "mass_case",
    "flux_case",
]

comp = const.merge(
    diss,
    on=merge_cols,
    suffixes=("_const", "_diss"),
)

comp["ratio_diss_const"] = (
    comp["weighted_mass_loss_kg_s_diss"] /
    comp["weighted_mass_loss_kg_s_const"]
)

comp["delta_log10_Mdot"] = np.log10(
    comp["ratio_diss_const"].replace(0, np.nan)
)

# Valid only if both lower-bound and upper-bound cases are still in Jeans regime
comp["both_jeans_valid"] = (
    comp["jeans_valid_const"] &
    comp["jeans_valid_diss"]
)

comp.to_csv("dissociation_vs_constant_comparison.csv", index=False)

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
    label=r"At least one case has $\lambda_J < 1.5$",
)

all_handles = atm_handles + flux_handles + [invalid_handle]


fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

for idx, mass in enumerate(mass_cases):
    ax = axes[idx]

    for atm in atm_archetype:
        for flux in flux_cases:

            s = comp[
                (comp["mass_case"] == mass) &
                (comp["atmosphere_type"] == atm) &
                (comp["flux_case"] == flux)
            ].sort_values("T_inf")

            if s.empty:
                continue

            valid = s[s["both_jeans_valid"]]
            invalid = s[~s["both_jeans_valid"]]

            ax.plot(
                valid["T_inf"],
                valid["delta_log10_Mdot"],
                color=species_colors[atm],
                linewidth=1.6,
                alpha=0.85,
            )

            ax.scatter(
                valid["T_inf"],
                valid["delta_log10_Mdot"],
                color=species_colors[atm],
                marker=flux_markers[flux],
                s=50,
            )

            ax.scatter(
                invalid["T_inf"],
                invalid["delta_log10_Mdot"],
                color=species_colors[atm],
                marker="x",
                s=85,
                linewidths=2,
            )

    ax.axhline(0, color="grey", linestyle="--", linewidth=1)

    ax.set_xlabel(r"$T_{\infty}$ [K]", fontsize=13)
    ax.set_title(mass_labels[mass], fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)

axes[0].set_ylabel(
    r"$\Delta \log_{10}\dot{M}_{\rm weighted}$ "
    r"$= \log_{10}(\dot{M}_{\rm diss}/\dot{M}_{\rm const})$",
    fontsize=12,
)

fig.legend(
    handles=all_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.12),
    ncol=4,
    frameon=False,
    fontsize=9,
)

plt.suptitle("Effect of full photodissociation on predicted Jeans escape", fontsize=15)
plt.tight_layout(rect=[0, 0.13, 1, 1])
plt.savefig(os.path.join(outdir, "dissociation_effect_ratio.png"), dpi=300, bbox_inches="tight")
plt.close()








#atmospheric lifetimes yaay

