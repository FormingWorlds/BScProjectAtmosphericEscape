import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.constants import G

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

case_styles = {
    "Constant VMR": "-",
    "Full dissociation": "--",
}


# Load lower and upper cases


constant = pd.read_csv("Outputs/proteus_jeans_case_summary.csv")
diss = pd.read_csv("Outputs/proteus_jeans_case_summary_dissociation.csv")

for table in [constant, diss]:
    if "unphysical_extension" not in table.columns:
        table["unphysical_extension"] = False

constant["case"] = "Constant VMR"
diss["case"] = "Full dissociation"

df = pd.concat([constant, diss], ignore_index=True)

df["log10_weighted_mass_loss"] = np.log10(
    df["weighted_mass_loss_kg_s"].replace(0, np.nan)
)


# Legend handles


atm_handles = [
    Line2D([0], [0], color=species_colors[a], lw=2, label=a)
    for a in atm_archetype
]

case_handles = [
    Line2D([0], [0], color="black", linestyle="-", label="Constant VMR"),
    Line2D([0], [0], color="black", linestyle="--", label="Full dissociation"),
]

flux_handles = [
    Line2D([0], [0], marker=flux_markers[f], color="black",
           linestyle="None", label=flux_labels[f])
    for f in flux_cases
]

hydro_handle = Line2D(
    [0], [0],
    marker="x",
    color="black",
    linestyle="None",
    markersize=8,
    markeredgewidth=2,
    label=r"Hydrodynamic onset ($\lambda_J < 1.5$)"
)


all_handles = atm_handles + case_handles + flux_handles + [
    hydro_handle]



def split_regimes(s):
    valid = s[
        (s["jeans_valid"]) &
        (~s["unphysical_extension"])
    ]

    hydro = s[
        (~s["jeans_valid"]) &
        (~s["unphysical_extension"])
    ]

    unphysical = s[s["unphysical_extension"]]

    return valid, hydro, unphysical



# 1. Escape bracket: constant VMR vs full dissociation


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

                valid, hydro, unphysical = split_regimes(s)

                ax.plot(
                    valid["T_inf"],
                    np.log10(valid["M_total_kg_s"]),
                    color=species_colors[atm],
                    linestyle=linestyle,
                    linewidth=1.7,
                    alpha=0.85,
                )

                ax.scatter(
                    valid["T_inf"],
                    np.log10(valid["M_total_kg_s"]),
                    color=species_colors[atm],
                    marker=flux_markers[flux],
                    s=50,
                )

                ax.scatter(
                    hydro["T_inf"],
                    np.log10(hydro["M_total_kg_s"]),
                    color=species_colors[atm],
                    marker="x",
                    s=80,
                    linewidths=2,
                )


    ax.set_xlabel(r"$T_{\infty}$ [K]", fontsize=13)
    ax.set_title(mass_labels[mass], fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)

axes[0].set_ylabel(
    r"$\log_{10}(\dot{M}_{\rm total})$ [kg/s]",
    fontsize=13,
)

fig.legend(
    handles=all_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.05),
    ncol=4,
    frameon=False,
    fontsize=9,
)

plt.tight_layout(rect=[0, 0.14, 1, 1])
plt.savefig(os.path.join(outdir, "escape_bracket_vs_Tinf.png"),
            dpi=300, bbox_inches="tight")
plt.close()



# 2. Dominant lambda transition plot


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

                valid, hydro, unphysical = split_regimes(s)

                physical = s[~s["unphysical_extension"]]

                ax.plot(
                    physical["T_inf"],
                    physical["dominant_lambda_j"],
                    color=species_colors[atm],
                    linestyle=linestyle,
                    linewidth=0.8,
                    alpha=0.55,
                )

                ax.scatter(
                    valid["T_inf"],
                    valid["dominant_lambda_j"],
                    color=species_colors[atm],
                    marker=flux_markers[flux],
                    s=35,
                )

                ax.scatter(
                    hydro["T_inf"],
                    hydro["dominant_lambda_j"],
                    color=species_colors[atm],
                    marker="x",
                    s=35,
                    linewidths=2,
                )


    ax.axhline(1.5, color="grey", linestyle="--", linewidth=1)
    ax.axhline(10, color="grey", linestyle=":", linewidth=1)
    #ax.set_ylim(1,10.5)

    ax.set_yscale("log")
    ax.set_xlabel(r"$T_{\infty}$ [K]", fontsize=13)
    ax.set_title(mass_labels[mass], fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)

axes[0].set_ylabel(
    r"Dominant species Jeans parameter $\lambda_J$",
    fontsize=13,
)

fig.legend(
    handles=all_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.05),
    ncol=4,
    frameon=False,
    fontsize=9,
)

#plt.suptitle("Transition from Jeans to non-Jeans escape", fontsize=15)
plt.tight_layout(rect=[0, 0.16, 1, 1])
plt.savefig(os.path.join(outdir, "dominant_lambda_transition_vs_Tinf.png"),
            dpi=300, bbox_inches="tight")
plt.close()


# 3. Lifetime plot

def read_bulk(comp, mass_case, flux_case):
    path = (
        f"PROTEUS data/{comp}_atmospheres/"
        f"planet_bulk_properties_{comp}_atmospheres_{mass_case}.csv"
    )

    bulk = pd.read_csv(path, sep="\t")
    row = bulk[bulk["Case"] == flux_case].iloc[0]

    return row["R_int [m]"], row["M_planet [kg]"]


P_SURF = 1e5
SECONDS_PER_YEAR = 365.25 * 24 * 3600


def one_bar_atmosphere_mass(R_p, M_p):
    g = G * M_p / R_p**2
    return 4 * np.pi * R_p**2 * P_SURF / g


M_atm_values = []

for _, row in df.iterrows():
    R_p, M_p = read_bulk(
        row["atmosphere_type"],
        row["mass_case"],
        row["flux_case"],
    )

    M_atm_values.append(one_bar_atmosphere_mass(R_p, M_p))

df["M_atm_1bar_kg"] = M_atm_values
df["lifetime_s"] = df["M_atm_1bar_kg"] / df["weighted_mass_loss_kg_s"]
df["lifetime_yr"] = df["lifetime_s"] / SECONDS_PER_YEAR
df["log10_lifetime_yr"] = np.log10(
    df["lifetime_yr"].replace([np.inf, 0], np.nan)
)

df.to_csv("Outputs/prediction_lifetimes_1bar.csv", index=False)

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

                valid, hydro, unphysical = split_regimes(s)

                ax.plot(
                    valid["T_inf"],
                    valid["log10_lifetime_yr"],
                    color=species_colors[atm],
                    linestyle=linestyle,
                    linewidth=1.5,
                    alpha=0.65,
                )

                ax.scatter(
                    valid["T_inf"],
                    valid["log10_lifetime_yr"],
                    color=species_colors[atm],
                    marker=flux_markers[flux],
                    s=35,
                )

                ax.scatter(
                    hydro["T_inf"],
                    hydro["log10_lifetime_yr"],
                    color=species_colors[atm],
                    marker="x",
                    s=35,
                    linewidths=2,
                )

    ax.axhline(6, color="grey", linestyle=":", linewidth=1)
    ax.axhline(9, color="grey", linestyle="--", linewidth=1)
    ax.axhline(np.log10(4.5e9), color="grey", linestyle="-.", linewidth=1)

    ax.text(7000, 6.1, "1 Myr", color="grey", fontsize=9, ha="right")
    ax.text(7000, 9.1, "1 Gyr", color="grey", fontsize=9, ha="right")
    ax.text(
        7000,
        np.log10(4.5e9) + 0.1,
        "4.5 Gyr",
        color="grey",
        fontsize=9,
        ha="right",
    )

    ax.set_xlabel(r"$T_{\infty}$ [K]", fontsize=13)
    ax.set_title(mass_labels[mass], fontsize=14)
    ax.tick_params(axis="both", which="major", labelsize=12)

axes[0].set_ylabel(
    r"$\log_{10}(\tau_{\rm loss,\,1bar})$ [yr]",
    fontsize=13,
)

# No unphysical handle here because they are not plotted in lifetime figure
lifetime_handles = atm_handles + case_handles + flux_handles + [hydro_handle]

fig.legend(
    handles=lifetime_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.05),
    ncol=4,
    frameon=False,
    fontsize=9,
)

#plt.suptitle("Predicted lifetime of a 1-bar equivalent atmosphere", fontsize=15)
plt.tight_layout(rect=[0, 0.14, 1, 1])
plt.savefig(os.path.join(outdir, "lifetime_1bar_vs_Tinf.png"),
            dpi=300, bbox_inches="tight")
plt.close()


# 4. Dissociation / constant VMR ratio plot


const = pd.read_csv("Outputs/proteus_jeans_case_summary.csv")
diss = pd.read_csv("Outputs/proteus_jeans_case_summary_dissociation.csv")

for table in [const, diss]:
    if "unphysical_extension" not in table.columns:
        table["unphysical_extension"] = False

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

comp["either_unphysical"] = (
    comp["unphysical_extension_const"] |
    comp["unphysical_extension_diss"]
)

comp["both_jeans_valid"] = (
    comp["jeans_valid_const"] &
    comp["jeans_valid_diss"] &
    (~comp["either_unphysical"])
)

comp.to_csv("Outputs/dissociation_vs_constant_comparison.csv", index=False)

ratio_handles = atm_handles + flux_handles + [
    Line2D(
        [0], [0],
        marker="x",
        color="black",
        linestyle="None",
        markersize=8,
        markeredgewidth=2,
        label=r"At least one case has $\lambda_J < 1.5$",
    ),
]

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

            hydro = s[
                (~s["both_jeans_valid"]) &
                (~s["either_unphysical"])
            ]

            unphysical = s[s["either_unphysical"]]

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
                hydro["T_inf"],
                hydro["delta_log10_Mdot"],
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
    handles=ratio_handles,
    loc="lower center",
    bbox_to_anchor=(0.5, -0.05),
    ncol=4,
    frameon=False,
    fontsize=9,
)

#plt.suptitle("Effect of full photodissociation on predicted Jeans escape", fontsize=15)
plt.tight_layout(rect=[0, 0.14, 1, 1])
plt.savefig(os.path.join(outdir, "dissociation_effect_ratio.png"),
            dpi=300, bbox_inches="tight")
plt.close()