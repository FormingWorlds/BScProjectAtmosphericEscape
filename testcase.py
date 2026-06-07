from jeans import *
from scipy.integrate import cumulative_trapezoid

# planet
R_earth = 6.371e6    # m
M_earth = 5.972e24   # kg


# species
m_O    = 16 * 1.6605e-27   # kg
masses = {"O": m_O}

r = np.linspace(R_earth, R_earth + 2e6, 500000)

# temperature — isothermal at exobase temperature 
T = np.full_like(r, 3000.0)

# hydrostatic density profile with varying g
g    = G * M_earth / r**2
H    = k * T / (m_O * g)
ln_n = cumulative_trapezoid(-1/H, r, initial=0)

n0 = 1e16

n = n0 * np.exp(ln_n)

species = {"O": n}

result = jeans_escape(
    r=r,
    T=T,
    species=species,
    species_masses=masses,
    M=M_earth
)


print("\n=== EXOBASE SUMMARY ===")
print(f"Altitude:    {result['exobase_altitude']/10**3:.2f} km")
print(f"Radius:      {result['exobase_radius']:.3e} m")
print(f"Temperature: {result['exobase_temperature']:.1f} K")

for sp, res in result["results"].items():
    print(f"\n--- {sp} ---")
    print(f"Lambda_J:                 {res['lambda_j']:.2f}")
    print(f"Mass loss rate:           {res['Mdot (kg/s)']:.3e} kg/s")
    print(f"Thermal velocity:         {res['v_th (m/s)']:.2f} m/s")
    print(f"Effusion velocity:        {res['effusion_velocity (m/s)']:.3e} m/s")
    print(f"Number density (exobase): {res['number_density (m^3)']:.3e} m⁻³")


m_H = 1 * 1.6605e-27

species = {
    "O": n,
    "H": 0.1 * n   
}

masses = {
    "O": m_O,
    "H": m_H
}

result = jeans_escape(
    r=r,
    T=T,
    species=species,
    species_masses=masses,
    M=M_earth
)

idx = result["exobase_index"]
r_e = result["exobase_radius"]
g_e = G * M_earth / r_e**2
m_mean_e = mean_mass(species, masses)[idx]
H_e = k * T[idx] / (m_mean_e * g_e)

sigma = 1e-19
n_tot = np.sum(list(species.values()), axis=0)

mfp_e = 1 / (sigma * n_tot[idx])   

print(f"altitude:  {(r_e - R_earth)/1e3:.2f} km")
print(f"g:         {g_e:.4f} m/s²")
print(f"H:         {H_e/1e3:.2f} km")
print(f"mfp:       {mfp_e/1e3:.2f} km")
print(f"mfp/H:     {mfp_e/H_e:.4f}  should be just above 1.0")

print("\n=== EXOBASE SUMMARY ===")
print(f"Altitude:    {result['exobase_altitude']/10**3:.2f} km")
print(f"Radius:      {result['exobase_radius']:.3e} m")
print(f"Temperature: {result['exobase_temperature']:.1f} K")

for sp, res in result["results"].items():
    print(f"\n--- {sp} ---")
    print(f"Lambda_J:                 {res['lambda_j']:.2f}")
    print(f"Mass loss rate:           {res['Mdot (kg/s)']:.3e} kg/s")
    print(f"Thermal velocity:         {res['v_th (m/s)']:.2f} m/s")
    print(f"Effusion velocity:        {res['effusion_velocity (m/s)']:.3e} m/s")
    print(f"Number density (exobase): {res['number_density (m^3)']:.3e} m⁻³")

