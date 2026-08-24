import numpy as np
import pandas as pd
from open_PROTEUS_csv import profiles, bulk_data
from open_PROTEUS_csv import elements, masses, fluxes


# Function to calculate fractional global mass loss:
# Has input options for isothermal/adiabatic atmosphere and v_imp/v_esc
# Returns X_loss vs m/M

def X_loss(v_V = 1, atm_type = "adiabatic"):
	"""Fractional mass loss of atmosphere
	   x = (v_imp*m) / (v_esc*M)       --> x-axis
	   * v_V = (v_imp)/(v_esc); assumes v_imp = v_esc as standard
	   * M = mass of planet impacted
	   * rho_pl = density of impactor
	   * atm_type: distinhuishes between atmosphere types:
	     - "isothermal"
	     - "adiabatic"
	"""
	m_M = np.linspace(0.0, 1.0, 100)
	x = v_V * m_M
	
	if atm_type == "isothermal":
	    X_loss = 0.4*x + 1.4*(x**2) - 0.8*(x**3)
	
	elif atm_type == "adiabatic":
	    X_loss = 0.4*x + 1.8*(x**2) - 1.2*(x**3)
	
	mask = (X_loss > 1)
	if mask.any():
	    i = np.argmax(mask)
	    X_loss = X_loss[:i]
	    m_M = m_M[:i]        
	
	return X_loss, m_M
    

# Function to convert m/M to radius for x-axis of global mass loss
    
def mM_to_r(m_M, M, rho_pl = 2000):
	"""Function to convert m/M to radius for giant impactor graphs
	   * M = mass of planet impacted
	   * rho_pl = density of impactor, standard as 2000 kg/m^3
	"""
	r_cubed = m_M * (3/4) * M / rho_pl
	r = (r_cubed)**(1/3)
	
	return r


# Getting fractional global mass loss for an adiabatic atmosphere and corresponding radius:

X_loss_adia, m_M = X_loss(v_V = 1, atm_type = "adiabatic")

# Adding global mass loss to the dictionary:

for e in elements:
    for m in masses:
        for f in fluxes:
            if e == "H2" and m == "1_M" and f == "1000_F":
                continue
            
            atm = profiles[e][m][f]
            bulk = bulk_data[e][m][f]
            
            M_loss_global = X_loss_adia * bulk["atm_mass"]
            r_imp_gl = mM_to_r(m_M, bulk["mass"], rho_pl = 2000)
            
            profiles[e][m][f]["M_loss_global"] = M_loss_global
            profiles[e][m][f]["m/M_global"] = m_M
            profiles[e][m][f]["r_imp_gl"] = r_imp_gl
            


            
# Adding the data to a csv file:

data = []

for e in elements:
    for m in masses:
        for f in fluxes:
            if e == "H2" and m == "1_M" and f == "1000_F":
                continue
            
            atm = profiles[e][m][f]
            bulk = bulk_data[e][m][f]

            filename = f"{e}_atmosphere_{m}_earth_{f}_earth.csv"

            r_imp_arr = atm["r_imp_gl"]
            Mloss_arr = atm["M_loss_global"]
            mM_arr = atm["m/M_global"]

            for rimp, Mloss, mM in zip(r_imp_arr, Mloss_arr, mM_arr):
                data.append({
                    "Filename": filename,
                    "Element": e,
                    "Planet mass": m,
                    "Earth Flux": f,
                    "Atmospheric mass [kg]": bulk["atm_mass"],
                    "Global mass loss [kg]": Mloss,
                    "m/M": mM,
                    "Radius impactor [m]": rimp
                })

df = pd.DataFrame(data)
df.to_csv("Outputs/Global_mass_loss_data.csv", index=False)




