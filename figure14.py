import numpy as np
from planetesimal_and_MLR import r_min, r_cap, r_gi
from open_PROTEUS_csv import profiles, bulk_data
from open_PROTEUS_csv import elements, masses, fluxes
import matplotlib.pyplot as plt

# Total impactor mass required:

def M_T(h, R, M_atm, M_planet, rho, rho_pl=2000):
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


for e in elements:
    for m in masses:
        for f in fluxes:
            if e == "H2" and m == "1_M" and f == "1000_F":
                continue

            atm  = profiles[e][m][f]
            bulk = bulk_data[e][m][f]

            # compute curve
            r, MT_Mp = M_T(
                bulk["h_avg"],
                bulk["radius"],
                bulk["atm_mass"],
                bulk["mass"],
                atm["rho"][0]
            )

            # get cutoffs
            rmin = bulk["r_min"]
            rcap = bulk["r_cap"]
            rgi  = bulk["r_gi"]

            # make figure
            plt.figure(figsize=(7,5))
            plt.style.use("tableau-colorblind10")

            plt.plot(r, MT_Mp, color='#006BA4', label=f"{e}, {m}, {f}")

            # vertical dashed lines
            plt.axvline(rmin, color='#FF800E', linestyle="--", alpha=0.7, label="r_min")
            plt.axvline(rcap, color='#ABABAB', linestyle="--", alpha=0.7, label="r_cap")
            plt.axvline(rgi, color='#595959', linestyle="--", alpha=0.7, label="r_gi")

            plt.xscale("log")
            plt.yscale("log")

            plt.xlabel("Impactor radius r [m]")
            plt.ylabel("Total impactor mass / planet mass (M_T / M_p)")
            plt.title(f"Total impactor mass needed to eject entire atmosphere ({e}, {m}, {f})")

            plt.grid(alpha=0.3)
            plt.legend()
            plt.tight_layout()
            
            save_name = f"MT_vs_r_{e}_{m}_{f}.png"
            plt.savefig(save_name, dpi=300, bbox_inches="tight")
            
            plt.show()
