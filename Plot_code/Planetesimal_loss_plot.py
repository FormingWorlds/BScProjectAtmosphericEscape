import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from open_PROTEUS_csv import elements, masses, fluxes

# Open data file
df = pd.read_csv("Outputs/Planetesimal_mass_loss_data.csv") 

# Colors for the elements
colors = {
    "H2": "#006BA4",
    "H2O": "#FF800E",
    "N2": "#ABABAB",
    "CO2": "#595959"
}

# Linestyle for fluxe
flux_style = {
    "1_F": "-",
    "1000_F": "--"
}

# Create figure:
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9,4), sharey=True)
plt.style.use("tableau-colorblind10")

# Left side: 1 M_earth
for element in elements:
    for flux in fluxes:
        subset = df[
            (df["Element"] == element) &
            (df["Planet mass"] == "1_M") &
            (df["Earth Flux"] == flux)
        ]

        if subset.empty:
            continue

        ax1.plot(
            subset["Radius impactor [m]"],
            subset["Planetesimal mass loss [kg]"],
            color=colors[element],
            linestyle=flux_style[flux]
        )

ax1.set_xlim(0, 1e6)
ax1.set_xlabel(r"Radius [m]")
ax1.set_ylabel(r"Mass ejected [kg]")
ax1.set_yscale("log")
ax1.set_title(r"$1\ M_{Earth}$")
ax1.grid(alpha=0.3)

# Right side: 10 M_Earth:
for element in elements:
    for flux in fluxes:
        subset = df[
            (df["Element"] == element) &
            (df["Planet mass"] == "10_M") &
            (df["Earth Flux"] == flux)
        ]

        if subset.empty:
            continue

        ax2.plot(
            subset["Radius impactor [m]"],
            subset["Planetesimal mass loss [kg]"],
            color=colors[element],
            linestyle=flux_style[flux]
        )

ax2.set_xlim(0, 2e7)
ax2.set_xlabel(r"Radius [m]")
ax2.set_title(r"$10\ M_{Earth}$")
ax2.grid(alpha=0.3)

# Shared legend:
flux_legend = [
    Line2D([0], [0], color="black", linestyle="-",  label=r"1 $F_{Earth}$"),
    Line2D([0], [0], color="black", linestyle="--", label=r"1000 $F_{Earth}$")
]

divider = Line2D([0], [0], color="white", linestyle="-", linewidth=1, label=" ")

elements_legend = [
    Line2D([0], [0], color="#006BA4", label="H2"),
    Line2D([0], [0], color="#FF800E", label="H2O"),
    Line2D([0], [0], color="#ABABAB", label="N2"),
    Line2D([0], [0], color="#595959", label="CO2")
]

fig.legend(
    handles = flux_legend + [divider] + elements_legend,
    loc="lower right",
    bbox_to_anchor=(1.15, 0.3),
    title="Legend"
)

fig.suptitle("Mass ejected by planetesimal impactors [kg] vs radii of impactors", y=1.03)
plt.tight_layout()
plt.savefig('Planetesimal mass loss.png', dpi=450, bbox_inches="tight")
plt.show()
