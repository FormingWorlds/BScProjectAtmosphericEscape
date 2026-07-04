import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from jeans import *
from scipy.integrate import cumulative_trapezoid
from scipy.constants import G, k
import numpy as np


#same test case as the oxygen isothermal atm
R_earth = 6.371e6
M_earth = 5.972e24

m_O = 16 * 1.6605e-27

masses = {
    "O": m_O
}

r = np.linspace(R_earth, R_earth + 2e6, 500000)

T = np.full_like(r, 3000.0)

g = G * M_earth / r**2
H = k * T / (m_O * g)

ln_n = cumulative_trapezoid(-1 / H, r, initial=0)

n0 = 1e16
n = n0 * np.exp(ln_n)

species = {
    "O": n
}


#dayside vs full
result_day = jeans_escape(
    r=r,
    T=T,
    species=species,
    species_masses=masses,
    M=M_earth,
    dayside=True
)

result_full = jeans_escape(
    r=r,
    T=T,
    species=species,
    species_masses=masses,
    M=M_earth,
    dayside=False
)

Mdot_day = result_day["results"]["O"]["Mdot (kg/s)"]
Mdot_full = result_full["results"]["O"]["Mdot (kg/s)"]

ratio = Mdot_full / Mdot_day


print("GEOMETRY TEST")

print(f"Dayside loss:      {Mdot_day:.6e} kg/s")
print(f"Full-surface loss: {Mdot_full:.6e} kg/s")

print(f"\nRatio full/day:    {ratio:.6f}")

if ratio == 2:
    print("Geometry test confirmed")
