from scipy.constants import G, k, pi, atomic_mass
import numpy as np
from scipy.integrate import cumulative_trapezoid
from jeans import *

# Constants

R_EARTH = 6.371e6 #m
M_EARTH = 5.972e24 #kg
RHO_EARTH = 5510 #kg/m^3


MASSES = {
    "N2": 28 * atomic_mass,
    "CO2": 44 * atomic_mass,
}

SIGMA = 1e-19  # fixed first, improve later with radius for each

FEUV_VALUES = [1, 2, 4, 6, 8, 10, 12, 14]
PLANET_MASSES = [0.8, 1.0, 1.2]

COMPOSITIONS = [
    {"CO2": 0.10, "N2": 0.90},
    {"CO2": 0.20, "N2": 0.80},
    {"CO2": 0.40, "N2": 0.60},
    {"CO2": 0.60, "N2": 0.40},
    {"CO2": 0.80, "N2": 0.20},
    {"CO2": 0.90, "N2": 0.10},
    {"CO2": 0.99, "N2": 0.01},
]


def radius_from_mass(M):
    return (3 * M / (4 * pi * RHO_EARTH)) ** (1 / 3)

def test_exobase_temperature(FEUV, X_CO2):
    T_lower = 1500.0

    # purely qualitative: EUV heats, CO2 cools
    heating = 15 * FEUV**2
    cooling = 50 * X_CO2

    return max(T_lower + heating - cooling, T_lower)


def test_atmosphere(Rp, Mp, T_exo, composition, n0=1.73e22, z_max=3e6, n_grid=60000):
    r = np.linspace(Rp, Rp + z_max, n_grid)
    T = np.full_like(r, T_exo)

    m_mean = sum(composition[sp] * MASSES[sp] for sp in composition)

    g = G * Mp / r**2
    H = k * T / (m_mean * g)

    ln_n = cumulative_trapezoid(-1 / H, r, initial=0)
    n_total = n0 * np.exp(ln_n)

    species = {
        sp: X * n_total
        for sp, X in composition.items()
    }

    return r, T, species


def weighted_mass_loss(result, composition):
    """
    Van Looveren et al. use an abundance-weighted average,
    not per species so added percentage for each element.
    """
    total = 0.0
    for sp, X in composition.items():
        total += X * result["results"][sp]["Mdot (kg/s)"]
    return total


rows = []

for M_factor in PLANET_MASSES:
    Mp = M_factor * M_EARTH
    Rp = radius_from_mass(Mp)

    for composition in COMPOSITIONS:
        X_CO2 = composition["CO2"]
        X_N2 = composition["N2"]

        for FEUV in FEUV_VALUES:
            T_exo = test_exobase_temperature(FEUV, X_CO2)
            #print(T_exo)

            try:
                r, T, species = test_atmosphere(
                    Rp=Rp,
                    Mp=Mp,
                    T_exo=T_exo,
                    composition=composition,
                )

                result = jeans_escape(
                    r=r,
                    T=T,
                    species=species,
                    species_masses=MASSES,
                    M=Mp,
                    sigma=SIGMA,
                    dayside=True,
                )

                mdot_weighted = weighted_mass_loss(result, composition)
                lambda_weighted = sum(
                    composition[sp] * result["results"][sp]["lambda_j"]
                    for sp in composition
                )

                rows.append({
                    "FEUV": FEUV,
                    "Mplanet_Mearth": M_factor,
                    "CO2_fraction": X_CO2,
                    "N2_fraction": X_N2,
                    "T_exo_K": result["exobase_temperature"],
                    "exobase_altitude_km": result["exobase_altitude"] / 1e3,
                    "Mdot_weighted_kg_s": mdot_weighted,
                    "log10_Mdot_weighted": np.log10(mdot_weighted) if mdot_weighted > 0 else -np.inf,
                    "lambda_weighted" : lambda_weighted,
                })

            except ValueError:
                rows.append({
                    "FEUV": FEUV,
                    "Mplanet_Mearth": M_factor,
                    "CO2_fraction": X_CO2,
                    "N2_fraction": X_N2,
                    "T_exo_K": T_exo,
                    "exobase_altitude_km": np.nan,
                    "Mdot_weighted_kg_s": np.nan,
                    "log10_Mdot_weighted": np.nan,
                })


df = pd.DataFrame(rows)

#print(df.head(20))

import matplotlib.pyplot as plt
#plt.style.use('science') 

composition_values = sorted(df["CO2_fraction"].unique())

cmap = plt.cm.Spectral_r            #to make it nice like in the paper 
colors = {
    x: cmap(i / (len(composition_values) - 1))
    for i, x in enumerate(composition_values)
}

#plotting T_exo

subset = df[df["Mplanet_Mearth"] == 1.0]

plt.figure()
for X_CO2 in sorted(subset["CO2_fraction"].unique()):
    s = subset[subset["CO2_fraction"] == X_CO2]
    plt.scatter(s["FEUV"], s["T_exo_K"], marker="o", linewidths=0.5, edgecolor ='black',
                color = colors[X_CO2], label=f"N2:{(1-X_CO2)*10**2:.0f}%, CO2:{X_CO2*10**2:.0f}%")

plt.xlabel(r"EUV flux [$F_{\mathrm{EUV},\oplus}$]")
plt.ylabel(r"Exobase temperature K")
plt.legend()
plt.grid(alpha=0.4)
plt.title("Trend validation: temperature exobase and EUV")
plt.tight_layout()
plt.savefig('Plots/exobase_EUV.png')
plt.show()

# Fig. 4-like: loss vs FEUV for 1 Earth mass
subset = df[df["Mplanet_Mearth"] == 1.0]

plt.figure()
for X_CO2 in sorted(subset["CO2_fraction"].unique()):
    s = subset[subset["CO2_fraction"] == X_CO2]
    plt.scatter(s["FEUV"], s["log10_Mdot_weighted"], marker="o", linewidths=0.5, edgecolor ='black',
                color = colors[X_CO2], label=f"N2:{(1-X_CO2)*10**2:.0f}%, CO2:{X_CO2*10**2:.0f}%")

plt.xlabel(r"EUV flux [$F_{\mathrm{EUV},\oplus}$]")
plt.ylabel(r"log10 weighted mass-loss rate [kg/s]")
plt.legend()
plt.grid(alpha=0.4)
plt.title("Trend validation: composition and EUV")
plt.tight_layout()
plt.savefig('Plots/trend_validation_mass_loss_EUV.png')
plt.show()

# Fig. 5-like: Jeans parameter vs FEUV for 1 Earth mass
plt.figure()
for X_CO2 in sorted(subset["CO2_fraction"].unique(), reverse = True):
    s = subset[subset["CO2_fraction"] == X_CO2]
    plt.scatter(s["FEUV"], np.log10(s["lambda_weighted"]), marker="o", linewidths=0.5, edgecolor ='black',
                color = colors[X_CO2], label=f"N2:{(1-X_CO2)*10**2:.0f}%, CO2:{X_CO2*10**2:.0f}%")

plt.axhline(np.log10(1.5), linestyle="--", label="blow-off limit")
plt.xlabel(r"EUV flux [$F_{\mathrm{EUV},\oplus}$]")
plt.ylabel(r"log10 Jeans parameter")
plt.grid(alpha = 0.4)
plt.legend()
plt.title("Trend validation: Jeans parameter")
plt.savefig('Plots/trend_validation_jeans_parameter.png')
plt.show()

# Fig. 6-like: mass variation for 99% CO2 / 1% N2
subset_mass = df[df["CO2_fraction"] == 0.99]

markers = ["v", "o", "^"]

plt.figure()
for i, M_factor in enumerate(sorted(subset_mass["Mplanet_Mearth"].unique())):
    s = subset_mass[subset_mass["Mplanet_Mearth"] == M_factor]
    plt.scatter(s["FEUV"], s["log10_Mdot_weighted"], marker=markers[i],linewidths=0.5, edgecolor ='black',
                label=f"{M_factor:.1f} Mearth")

plt.xlabel(r"EUV flux [$F_{\mathrm{EUV},\oplus}$]")
plt.ylabel(r"log10 weighted mass-loss rate [kg/s]")
plt.legend()
plt.grid(alpha = 0.4)
plt.title("Trend validation: planet mass")
plt.savefig('Plots/trend_validation_planet_mass.png')
plt.show()