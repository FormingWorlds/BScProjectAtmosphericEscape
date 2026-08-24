from figure14 import M_T
from planetesimal_and_MLR import r_min, r_cap, r_gi
import matplotlib.pyplot as plt

# Earth atmosphere and planet constants:

R_Earth = 6378e3	# Earth radius [m]
M_Earth = 5.98e24	# Earth mass [kg]
h_Earth = 8.5e3		# scaleheight [m]
M_atm_Earth = 5.15e18	# Earth atmosphere mass [kg]
rho_atm_Earth = 1.293	# density Earth atmosphere at 0 height [kg/m^3]
rho_pl = 2000		# density planetesimal [kg/m^3]

rmin = r_min(rho_atm_Earth, h_Earth, rho_pl=2000)
rcap = r_cap(rho_atm_Earth, h_Earth, R_Earth, rho_pl=2000)
rgi = r_gi(h_Earth, R_Earth)

# Plotting for Earth
r, MT_Mpl = M_T(h_Earth, R_Earth, M_atm_Earth, M_Earth, rho_atm_Earth, rho_pl=2000)

plt.figure(figsize=(6,5))
plt.style.use("tableau-colorblind10")

plt.plot(r, MT_Mpl)

plt.axvline(rmin, color='#FF800E', linestyle="--", alpha=0.7, label="r_min")
plt.axvline(rcap, color='#ABABAB', linestyle="--", alpha=0.7, label="r_cap")
plt.axvline(rgi, color='#595959', linestyle="--", alpha=0.7, label="r_gi")

plt.xscale("log")
plt.yscale("log")


plt.xlabel(r'$m/M$')
plt.ylabel(r'$X_{loss}$')
plt.title(r'Fractional global mass loss $X_{loss}$')
plt.grid(alpha=0.2)

plt.title("Total impactor mass needed to eject a full Earth-like atmosphere")

plt.legend(loc = 'lower right', bbox_to_anchor=(1.4, 0.1))

plt.savefig('Earth_total_M_imp_thresholds.png', dpi = 300, bbox_inches="tight")
plt.show()