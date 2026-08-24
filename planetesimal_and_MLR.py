import numpy as np
import pandas as pd
from open_PROTEUS_csv import profiles, bulk_data
from open_PROTEUS_csv import elements, masses, fluxes

# Impactor density:

rho_pl = 2000

# Define the threshold radii:

def r_min(rho, h, rho_pl=2000):
	"""Minimum radius of planetesimal impacts
	* rho = atmosphere density planet
	* rho_pl = density planetesimal impactor
	* h = scaleheight planet atmosphere
	"""
	r_min = (3*rho / rho_pl)**(1/3)*h
	return r_min


def r_cap(rho, h, R, rho_pl=2000):
	""" Cap size of planetesimal impactor radius
	* rho = atm. density planet
	* rho_pl = density planetesimal impactor
	* h = scaleheight planet atm.
	* R = radius planet
	"""
	r_cap = (3*np.sqrt(2*np.pi)*rho / (4*rho_pl))**(1/3) * (h*R)**(1/2)
	return r_cap


def r_gi(h, R):
	"""Radius where giant impacts dominate
	* h = scaleheight atm
	* R = radius Earth
	"""
	r_gi = (2*h*(R**2))**(1/3)
	return r_gi


# Calculate the maximum mass that can be ejected by a single planetesimal impactor:

def cap_mass(rho, h, R):
	"""Total cap mass for terrestial planets
	* rho = atmosphere density at the surface
	* h = (average) scale height of the atmosphere
	* R = radius of the planet
	"""
	M_cap = 2 * np.pi * rho * (h**2) * R
	return M_cap


# The ejected mass by a single planetesimal impactor:

def ejected_mass_planetesimal(h, R, rho, rho_pl = 2000):
	"""Ejected mass by a planetesimal impactor (r < r_gi)
	Once r > r_cap, the cap mass is the total ejected mass
	* h = scale height of atmosphere
	* R = radius of planet
	* rho = atmosphere density at surface
	* rho_pl = impactor density; set to 2000 kg/m^3
	"""
	r_min_val = r_min(rho, h, rho_pl)
	r_cap_val = r_cap(rho, h, R, rho_pl)
	r_gi_val = r_gi(h, R)

	r_range = np.linspace(r_min_val, r_gi_val, 1000)
	M_ejected = []

	for r in r_range:

		if r < r_min_val:
			M_eject = 0.0
			print("No mass loss possible due to too small impactor size.")

		elif r >= r_min_val and r < r_cap_val:
			m_imp = rho_pl * (4/3) * np.pi * r**3
			M_eject = (r_min_val/(2*r)) * (1 - (r_min_val/r)**2) * m_imp

		elif r >= r_cap_val:
			M_eject = cap_mass(rho, h, R)

		M_ejected.append(M_eject)

	M_ejected = np.array(M_ejected)

	return M_ejected, r_range


# Calculate the mass loss rate in kg/s

def mass_loss_rate(q, M_pl, h, R, rho, rho_pl = 2000):
	"""Mass loss rate dM_atm/dt in kg/s
    * q = differential power law index
    * M_pl = total impactor mass per unit time
    * h = scale height
    * R = radius of planet
    * rho = atmosphere density at surface
    * rho_pl = density of impactor; set to 2000 kg/m^3
    """
	# Ejected mass by planetesimal impactors
	M_eject, r_range = ejected_mass_planetesimal(h, R, rho, rho_pl)

	# Mass of the impactor:
	m_pl = rho_pl * (4/3) * np.pi * (r_range)**3

	mass_int = (r_range)**(-q) * M_eject
	Mass_I = np.trapz(mass_int, r_range)

	m_int = (r_range)**(-q) * m_pl
	m_I = np.trapz(m_int, r_range)
    
	ratio = Mass_I / m_I
	# print("Integral ratio =", ratio)

	# Mass loss rate in kg/s (without minus sign: positive mass loss)
	dM_dt = M_pl * ratio

	return dM_dt


# Adding threshold radii, cap mass, and planetesimal mass loss to the dictionary:

for e in elements:
    for m in masses:
        for f in fluxes:
            if e == "H2" and m == "1_M" and f == "1000_F":
                continue
            
            atm = profiles[e][m][f]
            bulk = bulk_data[e][m][f]
            
            r_min_val = r_min(atm["rho"][0], bulk["h_avg"], rho_pl = 2000)
            r_cap_val = r_cap(atm["rho"][0], bulk["h_avg"], bulk["radius"], rho_pl=2000)
            r_gi_val = r_gi(bulk["h_avg"], bulk["radius"])
            
            bulk_data[e][m][f]["r_min"] = r_min_val
            bulk_data[e][m][f]["r_cap"] = r_cap_val
            bulk_data[e][m][f]["r_gi"] = r_gi_val
            
            M_cap = cap_mass(atm["rho"][0], bulk["h_avg"], bulk["radius"])
            bulk_data[e][m][f]["M_cap"] = M_cap
            
            M_eject, r = ejected_mass_planetesimal(bulk["h_avg"], bulk["radius"], atm["rho"][0], rho_pl = 2000)
            profiles[e][m][f]["M_eject_pl"] = M_eject
            profiles[e][m][f]["r_imp_pl"] = r
            
            
# Adding the data for planetesimal mass loss to a csv file:

data = []

for e in elements:
    for m in masses:
        for f in fluxes:
            if e == "H2" and m == "1_M" and f == "1000_F":
                continue
            
            atm = profiles[e][m][f]
            bulk = bulk_data[e][m][f]

            filename = f"{e}_atmosphere_{m}_earth_{f}_earth.csv"

            r_imp_arr = atm["r_imp_pl"]
            Mloss_arr = atm["M_eject_pl"]

            for rimp, Mloss in zip(r_imp_arr, Mloss_arr):
                data.append({
                    "Filename": filename,
                    "Element": e,
                    "Planet mass": m,
                    "Earth Flux": f,
                    "r_min [m]": bulk["r_min"],
                    "r_cap [m]": bulk["r_cap"],
                    "r_gi [m]": bulk["r_gi"],
                    "rho_0 [kg/m^3]": atm["rho"][0],
                    "Cap mass [kg]": bulk["M_cap"],
                    "Atmospheric mass [kg]": bulk["atm_mass"],
                    "Radius impactor [m]": rimp,
                    "Planetesimal mass loss [kg]": Mloss
                })

df = pd.DataFrame(data)
df.to_csv("Outputs/Planetesimal_mass_loss_data.csv", index=False)


# Defining ranges for mass loss rate calculations:
# Set q, different M_pl:
q_set= 3.0
M_pl_range = np.logspace(2, 13, 500)

# Different q, set M_pl:
q_range = np.linspace(1.1, 4.0, 500)
M_pl_set = 10**7


# Calculating mass loss rate vs. total impactor mass per time:
for e in elements:
    for m in masses:
        for f in fluxes:
            if e == "H2" and m == "1_M" and f == "1000_F":
                continue
            
            atm = profiles[e][m][f]
            bulk = bulk_data[e][m][f]
            
            MLR_vs_Mpl = mass_loss_rate(q_set, M_pl_range, bulk["h_avg"], bulk["radius"], atm["rho"][0], rho_pl = 2000)
            MLR_vs_q = [mass_loss_rate(q, M_pl_set, bulk["h_avg"], bulk["radius"], atm["rho"][0], rho_pl = 2000) for q in q_range]
            
            # dM/dt for set q and changing impactor mass
            profiles[e][m][f]["MLR_per_Mpl"] = MLR_vs_Mpl
            profiles[e][m][f]["M_pl_range"] = M_pl_range
            bulk_data[e][m][f]["Set_q"] = q_set
            
            # dM/dt for varying q and set M_pl:
            profiles[e][m][f]["MLR_per_q"] = MLR_vs_q
            profiles[e][m][f]["q_range"] = q_range
            bulk_data[e][m][f]["Set_M_pl"] = M_pl_set


# Adding the data for MLR per M_pl to a csv file:

data2 = []

for e in elements:
    for m in masses:
        for f in fluxes:
            if e == "H2" and m == "1_M" and f == "1000_F":
                continue
            
            atm = profiles[e][m][f]
            bulk = bulk_data[e][m][f]

            filename = f"{e}_atmosphere_{m}_earth_{f}_earth.csv"

            Mpl_arr = atm["M_pl_range"]
            MLR_arr = atm["MLR_per_Mpl"]

            for Mpl, MLR in zip(Mpl_arr, MLR_arr):
                data2.append({
                    "Filename": filename,
                    "Element": e,
                    "Planet mass": m,
                    "Earth Flux": f,
                    "r_min [m]": bulk["r_min"],
                    "r_cap [m]": bulk["r_cap"],
                    "r_gi [m]": bulk["r_gi"],
                    "Cap mass [kg]": bulk["M_cap"],
                    "Differential power law index q": bulk["Set_q"],
                    "Total impactor mass [kg/s]": Mpl,
                    "Mass loss rate [kg/s]": MLR
                })

df = pd.DataFrame(data2)
df.to_csv("Outputs/Mass_loss_rate_vs_M_pl.csv", index=False)


# Adding the data for MLR per q to a csv file:

data3 = []

for e in elements:
    for m in masses:
        for f in fluxes:
            if e == "H2" and m == "1_M" and f == "1000_F":
                continue
            
            atm = profiles[e][m][f]
            bulk = bulk_data[e][m][f]

            filename = f"{e}_atmosphere_{m}_earth_{f}_earth.csv"

            q_arr = atm["q_range"]
            MLR_arr = atm["MLR_per_q"]

            for q, MLR in zip(q_arr, MLR_arr):
                data3.append({
                    "Filename": filename,
                    "Element": e,
                    "Planet mass": m,
                    "Earth Flux": f,
                    "r_min [m]": bulk["r_min"],
                    "r_cap [m]": bulk["r_cap"],
                    "r_gi [m]": bulk["r_gi"],
                    "Cap mass [kg]": bulk["M_cap"],
                    "Total impactor mass [kg]/s": bulk["Set_M_pl"],
                    "Differential power law index q": q,
                    "Mass loss rate [kg/s]": MLR,
                })

df = pd.DataFrame(data3)
df.to_csv("Outputs/Mass_loss_rate_vs_q.csv", index=False)

