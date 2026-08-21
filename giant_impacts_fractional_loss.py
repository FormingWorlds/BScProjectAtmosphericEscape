import numpy as np


def X_loss(gamma, v_V = 1):
        """Fractional mass loss of atmosphere
        x = (v_imp*m) / (v_esc*M)       --> x-axis
        * v_V = (v_imp)/(v_esc); assumes v_imp = v_esc as standard
        * Gamma: adiabatic index
                 - 5/3 (isothermal)
                 - 4/3 (adiabatic)
        * M = mass of planet impacted
        * rho_pl = density of impactor
        """
        m_M = np.linspace(0.0, 1.0, 100)
        x = v_V * m_M
	
        if gamma == 5/3:
                X_loss = 0.4*x + 1.4*(x**2) - 0.8*(x**3)

        elif gamma == 4/3:
                X_loss = 0.4*x + 1.8*(x**2) - 1.2*(x**3)

        mask = (X_loss > 1)
        if mask.any():
            i = np.argmax(mask)
            X_loss = X_loss[:i]
            m_M = m_M[:i]
            
        
        return X_loss, m_M
    

    
def mM_to_r(m_M, M, rho_pl = 2000):
    """Function to convert m/M to radius for giant impactor graphs
       * M = mass of planet impacted
       * rho_pl = density of impactor, standard as 2000 kg/m^3
    """
    r_cubed = m_M * (3/4) * M / rho_pl
    r = (r_cubed)**(1/3)
    
    return r

