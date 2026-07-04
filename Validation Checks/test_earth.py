import numpy as np
from scipy.constants import G, k, pi, atomic_mass
import pandas as pd

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from jeans import jeans_escape

R_EARTH = 6.371e6       # m
M_EARTH = 5.972e24      # kg

SPECIES_MASSES = {
    "H": 1.0 * atomic_mass,
}

# Approximate present-day Earth exobase conditions
z_exo = 500e3           # m, exobase altitude
T_exo = 1000.0          # K, thermospheric/exospheric temperature

# For Earth, eq. (5.13) gives ve ~10.8 km s–1 at the exobase (Catling_and_Kasting_ch_5). 
# Thus, nexob 3.6e07 cm–3, which is roughly similar to the value at the exobase on most planets and satellites.(including oxygen so should use smaller number)
#The number density at Earth’s exobase can be estimated from eq. (5.9), taking σc ~3*10^-15 cm2

# H density chosen for present-day Earth order-of-magnitude Jeans escape.

n_H_exo = 1e11        # m^-3 where to get this from?

sigma = 3e-19           # m^2, approximate


def make_earth_h_profile():
    """
    Construct a minimal hydrogen exosphere profile.

    This is not a full Earth atmosphere model, only includes hydrogen.
    It starts at the assumed exobase and tests Jeans escape mass loss.
    """
    r = np.linspace(R_EARTH + z_exo, R_EARTH + 100e5, 100)
    T = np.full_like(r, T_exo)

    species = {
        "H": np.full_like(r, n_H_exo)
    }

    return r, T, species


r, T, species = make_earth_h_profile()

result = jeans_escape(
    r=r,
    T=T,
    species=species,
    species_masses=SPECIES_MASSES,
    M=M_EARTH,
    sigma=sigma,
    dayside=False,   # global Earth area, 4*pi*r^2
)

H = result["results"]["H"]

print("=== Earth hydrogen Jeans escape test ===")
print(f"Exobase altitude: {result['exobase_altitude'] / 1e3:.1f} km")
print(f"Exobase temperature: {result['exobase_temperature']:.1f} K")
print(f"Hydrogen density at exobase: {H['number_density (m^3)']:.3e} m^-3")
print(f"lambda_J(H): {H['lambda_j']:.3f}")
print(f"v_th(H): {H['v_th (m/s)']:.3e} m/s")
print(f"effusion velocity: {H['effusion_velocity (m/s)']:.3e} m/s")
print(f"Mdot_H: {H['Mdot (kg/s)']:.3e} kg/s")

target = 3 #kg/s
factor = H["Mdot (kg/s)"] / target

print(f"Target scale: ~3 kg/s")
print(f"Model / target: {factor:.2f}")

if 0.1 < factor < 10:
    print("Confirmed: Earth H escape is within an order of magnitude.")
else:
    print("Re-check: Earth H escape is outside the expected scale.")

physical_altitude = (result["exobase_radius"] - R_EARTH)/1e3 #this is because my grid starts at exobase so it will show 0 km

print(
    f"Physical exobase altitude: "
    f"{physical_altitude:.1f} km"
)

n_H_values = [1e6, 3e6, 1e7, 3e7, 1e8, 1e9, 1e10, 1e11]  # m^-3
T_values = [700, 900, 1000, 1200, 1500]

rows = []

for T_exo in T_values:
    r = np.linspace(R_EARTH + z_exo, R_EARTH + z_exo + 100e3, 100)
    T = np.full_like(r, T_exo)

    for n_H_exo in n_H_values:
        species = {
            "H": np.full_like(r, n_H_exo)
        }

        result = jeans_escape(
            r=r,
            T=T,
            species=species,
            species_masses=SPECIES_MASSES,
            M=M_EARTH,
            sigma=sigma,
            dayside=False,
        )

        H = result["results"]["H"]

        rows.append({
            "T_exo_K": T_exo,
            "n_H_m3": n_H_exo,
            "lambda_J": H["lambda_j"],
            "Mdot_H_kg_s": H["Mdot (kg/s)"],
            "model_over_target": H["Mdot (kg/s)"] / 1e-4,
        })

df = pd.DataFrame(rows)
print(df)

import matplotlib.pyplot as plt
import numpy as np

CB_color_cycle = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

plt.figure(figsize=(7, 5))

for i, T_exo in enumerate(T_values):
    s = df[df["T_exo_K"] == T_exo]
    plt.plot(
        s["n_H_m3"],
        s["Mdot_H_kg_s"],
        marker="o",
        label=f"{T_exo} K",
        color=CB_color_cycle[i]
        )

plt.axhline(3, color="black", linestyle="--", label=r"$3$ kg/s target")

plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$n_{\rm H,exo}$ [m$^{-3}$]")
plt.ylabel(r"$\dot{M}_{\rm H}$ [kg/s]")
plt.title("Earth H Jeans escape validation sensitivity")
plt.legend()
plt.tight_layout()
plt.savefig("Plots/earth_H_validation_sensitivity.png", dpi=300, bbox_inches="tight")
plt.show()