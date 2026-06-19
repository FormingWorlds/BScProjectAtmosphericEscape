import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

outdir = "Plots/Dissociation_comparison"
os.makedirs(outdir, exist_ok=True)

T_choice = 1000
mass_choice = "10_M_earth"
flux_choice = "1000_F_earth"

atm_order = ["CO2", "N2", "H2O", "H2"]

atm_labels = {
    "CO2": r"CO$_2$",
    "N2": r"N$_2$",
    "H2O": r"H$_2$O",
    "H2": r"H$_2$",
}

constant = pd.read_csv("proteus_jeans_case_summary.csv")
diss = pd.read_csv("proteus_jeans_case_summary_dissociation.csv")

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

# =========================
# Option A: grouped bar chart
# =========================

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
    color="#56B4E9"
)

for i, atm in enumerate(atm_order):

    const_val = pivot.loc[atm, "Constant VMR"]
    diss_val = pivot.loc[atm, "Full dissociation"]

    const_species = s[
        (s["atmosphere_type"] == atm) &
        (s["case"] == "Constant VMR")
    ]["dominant_escaping_species"].iloc[0]

    diss_species = s[
        (s["atmosphere_type"] == atm) &
        (s["case"] == "Full dissociation")
    ]["dominant_escaping_species"].iloc[0]

    ax.text(
        x[i] - width / 2,
        np.log10(const_val) + 0.05,
        const_species,
        ha="center",
        va="bottom",
        fontsize=9,
        rotation=0,
    )

    ax.text(
        x[i] + width / 2,
        np.log10(diss_val) + 0.05,
        diss_species,
        ha="center",
        va="bottom",
        fontsize=9,
        rotation=0,
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


# =========================
# Option B: scatter comparison
# =========================

comp = constant.merge(
    diss,
    on=["file", "T_inf", "atmosphere_type", "mass_case", "flux_case"],
    suffixes=("_const", "_diss")
)

comp = comp[
    (comp["T_inf"] == T_choice) &
    (comp["mass_case"] == mass_choice) &
    (comp["flux_case"] == flux_choice) &
    (~comp["unphysical_extension_const"]) &
    (~comp["unphysical_extension_diss"])
].copy()

fig, ax = plt.subplots(figsize=(6, 6))

for atm in atm_order:
    a = comp[comp["atmosphere_type"] == atm]

    if a.empty:
        continue

    xval = a["dominant_species_Mdot_kg_s_const"].iloc[0]
    yval = a["dominant_species_Mdot_kg_s_diss"].iloc[0]

    ax.scatter(
        np.log10(xval),
        np.log10(yval),
        s=80,
        label=atm_labels[atm]
    )

    ax.annotate(
        atm_labels[atm],
        (np.log10(xval), np.log10(yval)),
        xytext=(5, 5),
        textcoords="offset points",
        fontsize=10
    )

lims = [
    min(
        np.log10(comp["dominant_species_Mdot_kg_s_const"]).min(),
        np.log10(comp["dominant_species_Mdot_kg_s_diss"]).min()
    ) - 0.5,
    max(
        np.log10(comp["dominant_species_Mdot_kg_s_const"]).max(),
        np.log10(comp["dominant_species_Mdot_kg_s_diss"]).max()
    ) + 0.5
]

ax.plot(
    lims,
    lims,
    color="grey",
    linestyle="--",
    linewidth=1,
    label="No change"
)

ax.set_xlim(lims)
ax.set_ylim(lims)

ax.set_xlabel(
    r"$\log_{10}(\dot{M}_{\rm dominant, const})$ [kg s$^{-1}$]"
)
ax.set_ylabel(
    r"$\log_{10}(\dot{M}_{\rm dominant, diss})$ [kg s$^{-1}$]"
)

ax.set_title(
    fr"Effect of full dissociation at $T_\infty={T_choice}$ K"
)

ax.legend(frameon=False)
plt.tight_layout()

plt.savefig(
    os.path.join(outdir, "dominant_Mdot_scatter_constant_vs_dissociation.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.close()