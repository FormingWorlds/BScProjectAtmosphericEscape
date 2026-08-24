import numpy as np
from planetesimal_and_MLR import r_min, r_cap, r_gi
from M_imp_tot_per_rad import M_T_plot
from open_PROTEUS_csv import profiles, bulk_data
from open_PROTEUS_csv import elements, masses, fluxes
import matplotlib.pyplot as plt

# Add the data of M_imp per radius to the dictionary for the models

for e in elements:
    for m in masses:
        for f in fluxes:
            if e == "H2" and m == "1_M" and f == "1000_F":
                continue

            atm  = profiles[e][m][f]
            bulk = bulk_data[e][m][f]

            # compute curve
            r, MT_Mp = M_T_plot(
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
