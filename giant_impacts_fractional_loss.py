import numpy as np

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
    
    
# Options for mass of planet:

M_1Earth = 5.92e24
M_10Earth = 10*M_1Earth
    

# Function to convert m/M to radius for x-axis of global mass loss
    
def mM_to_r(m_M, M, rho_pl = 2000):
	"""Function to convert m/M to radius for giant impactor graphs
	   * M = mass of planet impacted
	   * rho_pl = density of impactor, standard as 2000 kg/m^3
	"""
	r_cubed = m_M * (3/4) * M / rho_pl
	r = (r_cubed)**(1/3)
	
	return r

