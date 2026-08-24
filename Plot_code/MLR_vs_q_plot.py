import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from open_PROTEUS_csv import elements, masses, fluxes

# Get data from csv file
df = pd.read_csv("Outputs/Mass_loss_rate_vs_q.csv")   # <-- your file name

# Colors for each element
colors = {
    "H2": "#006BA4",
    "H2O": "#FF800E",
    "N2": "#ABABAB",
    "CO2": "#595959"
}

# Linestyle for fluxes
flux_style = {
    "1_F": "-",
    "1000_F": "--"
}

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4.5), sharey=True)
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
            subset["Differential power law index q"],
            subset["Mass loss rate [kg/s]"],
            color=colors[element],
            linestyle=flux_style[flux]
        )

ax1.set_yscale("log")
ax1.set_xlabel(r"Differential power law index q")
ax1.set_ylabel(r"|$dM_{atm}/dt$| [kg/s]")
ax1.set_title("1 M$_{Earth}$")
ax1.grid(alpha=0.3)

# Right side: 10 M_earth
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
            subset["Differential power law index q"],
            subset["Mass loss rate [kg/s]"],
            color=colors[element],
            linestyle=flux_style[flux]
        )

ax2.set_yscale("log")
ax2.set_xlabel(r"Differential power law index q")
ax2.set_title("10 M$_{Earth}$")
ax2.grid(alpha=0.3)

# Legend
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
    loc="center right",
    bbox_to_anchor=(1.17, 0.5),
    title="Legend"
)

fig.suptitle(r"Mass loss rate vs differential power law index q for $M_{pl}=10^7\ \mathrm{kg/s}$")
plt.tight_layout()
plt.savefig('MLR_vs_q.png', dpi=450, bbox_inches="tight")
plt.show()
