import numpy as np
import matplotlib.pyplot as plt

r_min = 1       # minimum impact size (km)
r_cap = 25      # cap radius (km)

gamma_iso = 5/3 # adiabatic index for isothermal atmosphere
beta_iso = 1.90

rho_0 = 1.293   # density on the ground (kg/m^3)
h_earth = 8.5   # scale height earth (km)



m = 6.4e23      # mass of impactor: Mars (kg)
M = 5.98e24     # mass of planet: Earth (kg)


# v_ratio = (v_imp/v_esc)*(m/M)
v_ratio = np.linspace(0.0, 1.0, 100)
v_equ = m/M
print(v_equ)

X_loss_iso = 0.4*(v_ratio) + 1.4*(v_ratio)**2 - 0.8*(v_ratio)**3


plt.plot(v_ratio, X_loss_iso, color='red')
plt.axvline(v_equ, color='navy')
plt.xlabel('(v_imp/v_esc)(m/M)')
plt.ylabel('X_loss')
plt.title('X_loss from Mars-like giant impact on isothermal Earth')
#plt.xscale('log')
#plt.yscale('log')
plt.show()
#plt.savefig(isothermal_model_giant)
