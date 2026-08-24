import numpy as np
from planetesimal_and_MLR import r_min, r_cap, r_gi

# Total impactor mass required:

def M_T_plot(h, R, M_atm, M_planet, rho, rho_pl=2000):
    """Total impactor mass required for different impactor radii, based on Schlichting et al. (2015)
    * h = scale height [m]
    * R = radius planet [m]
    * M_atm = atmosphere mass [kg]
    * M_planet = planet mass [kg]
    * rho = surface density [kg/m^3]
    """
    
    r_min_val = r_min(rho, h, rho_pl)
    r_cap_val = r_cap(rho, h, R, rho_pl)
    r_gi_val  = r_gi(h, R)

    # radius range (log-spaced)
    r = np.geomspace(r_min_val + 1e-6, 3e8, 2000)

    # initialize output
    M_T = np.zeros_like(r)

    # region 1: r < r_cap
    mask1 = r < r_cap_val
    M_T[mask1] = (2*r[mask1] / r_min_val) * (1 - (r_min_val/r[mask1])**2)**(-1) * M_atm

    # region 2: r_cap <= r < r_gi
    mask2 = (r >= r_cap_val) & (r < r_gi_val)
    M_T[mask2] = (4*np.pi/3) * rho_pl * r[mask2]**3 * (2*R/h)

    # region 3: r >= r_gi
    mask3 = r >= r_gi_val
    M_T[mask3] = 4 * M_planet

    # normalize
    MT_Mplanet = M_T / M_planet
    return r, MT_Mplanet