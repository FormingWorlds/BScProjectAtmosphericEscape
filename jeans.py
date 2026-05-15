from physics import *
from exobase import *
from cross_section import *

# escape rate
def jeans_escape(r, T, species, species_masses, M, sigma = 1e-19, dayside=True):
    
    n_tot = sum(species.values())
    m_mean = mean_mass(species, species_masses)

    if sigma == "weighted":
        sigma_used = effective_cross_section(species)
    else:
        sigma_used = sigma
    
    idx = find_exobase(r, n_tot, T, m_mean, M, sigma_used)

    sigma_exo = sigma_used[idx] if np.ndim(sigma_used) > 0 else sigma_used
    
    r_exo, T_exo = r[idx], T[idx]
    A = (2 if dayside else 4) * pi * r_exo**2

    exobase_altitude = r_exo - r[0]
    
    results = {}

    M_weighted = 0
    
    for sp, n in species.items():
        m = species_masses[sp]
        n_exo = n[idx]
        
        v_th = thermal_velocity(T_exo, m)
        lam = jeans_parameter(M, m, T_exo, r_exo)
        f_j = effusion_velocity(v_th, lam)
        
        Mdot = n_exo * f_j * A * m
        X_exo = n_exo / n_tot[idx]
        M_weighted += X_exo * Mdot
        
        results[sp] = {
            "Mdot (kg/s)": float(Mdot),
            "lambda_j": float(lam),
            "v_th (m/s)": float(v_th),
            "effusion_velocity (m/s)": float(f_j),
            "number_density (m^3)": float(n_exo)
        }
    
    return {
        "weighted_mass_loss_kg_s": M_weighted,
        "results": results,
        "exobase_index": int(idx) , 
        "exobase_radius": float(r_exo),
        "temperature": float(T_exo), 
        "exobase_altitude": float(exobase_altitude),
        "sigma_exobase_m2": float(sigma_exo),
    }