import numpy as np

# x = (v_imp*m) / (v_esc*M) --> x-axis of graphs


def X_loss(gamma, v_V):
        """Fractional mass loss of atmosphere
        x = (v_imp*m) / (v_esc*M)       --> x-axis
        v_V = (v_imp)/(v_esc)
        Gamma: adiabatic index
                - 5/3 (isothermal)
                - 4/3 (adiabatic)
        """
        m_M = np.linspace(0.0, 1.0, 100)
        x = v_V * m_M
	
        if gamma == 5/3:
                X_loss = 0.4*x + 1.4*(x**2) - 0.8*(x**3)

        elif gamma == 4/3:
                X_loss = 0.4*x + 1.8*(x**2) - 1.2*(x**3)

        return X_loss, m_M


