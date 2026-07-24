import numpy as np

# x = (v_imp*m) / (v_esc*M) --> x-axis of graphs

x = np.linspace(0.0, 1.0, 100)

def X_loss(x, gamma):
        """Fractional mass loss of atmosphere
        x = (v_imp*m) / (v_esc*M)       --> x-axis
        Gamma: adiabatic index
                - 5/3 (isothermal)
                - 4/3 (adiabatic)
        """
        if gamma == 5/3:
                X_loss = 0.4*x + 1.4*(x**2) - 0.8*(x**3)

        elif gamma == 4/3:
                X_loss = 0.4*x + 1.8*(x**2) - 1.2*(x**3)

        return X_loss


