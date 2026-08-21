import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from giant_impacts_fractional_loss import X_loss


gamma_adia = 4/3 # index: adiabatic
gamma_iso = 5/3 # index: isothermal


# values of v_imp/v_esc used in Schlichting:
v_V = [3, 2, 1.5, 1.2, 1, 0.8, 0.5, 0.2]

# colors for different values of v_imp/v_esc:
colors = ['#006BA4', '#FF800E', '#ABABAB', '#595959', '#5F9ED1', '#C85200', '#898989', '#A2C8EC']


# Start the plot and X_loss calculation:

plt.figure(figsize=(6,5))
plt.style.use("tableau-colorblind10")

for v, c in zip(v_V, colors):
    X_loss_is, m_M_is = X_loss(gamma_iso, v) # Isothermal mass loss
    X_loss_ad, m_M_ad = X_loss(gamma_adia, v) # Adiabatic mass loss
    
    plt.plot(m_M_ad, X_loss_ad, color = c)
    plt.plot(m_M_is, X_loss_is, color = c, linestyle='--')

plt.xlabel(r'$m/M$')
plt.ylabel(r'$X_{loss}$')
plt.title(r'Fractional global mass loss $X_{loss}$')
plt.grid(alpha=0.2)

# Legend explaining the different colors and dashed vs. solid lines
linestyle_legend = [
    Line2D([0], [0], color="black", linestyle="-",  label="Adiabatic"),
    Line2D([0], [0], color="black", linestyle="--", label="Isothermal")
]

divider = Line2D([0], [0], color="white", linestyle="-", linewidth=1, label=" ")

velocity_legend = [
    Line2D([0], [0], color="#006BA4", label=r"$v_{imp}/v_{esc}=$3"),
    Line2D([0], [0], color="#FF800E",   label=r"$v_{imp}/v_{esc}=$2"),
    Line2D([0], [0], color="#ABABAB",  label=r"$v_{imp}/v_{esc}=$1.5"),
    Line2D([0], [0], color="#595959",  label=r"$v_{imp}/v_{esc}=$1.2"),
    Line2D([0], [0], color="#5F9ED1",  label=r"$v_{imp}/v_{esc}=$1"),
    Line2D([0], [0], color="#C85200",  label=r"$v_{imp}/v_{esc}=$0.8"),
    Line2D([0], [0], color="#898989",  label=r"$v_{imp}/v_{esc}=$0.5"),
    Line2D([0], [0], color="#A2C8EC", label=r"$v_{imp}/v_{esc}=$0.2")
]

plt.legend(handles = linestyle_legend + [divider] + velocity_legend, loc = 'lower right', dpi = 300, bbox_to_anchor=(1.4, 0.1))

plt.savefig('Fractional_global_loss.png', bbox_inches="tight")
plt.show()


