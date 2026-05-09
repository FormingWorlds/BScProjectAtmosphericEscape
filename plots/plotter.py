#functions for plotting parameters of the planets and the escape stuff
import os
import matplotlib.pyplot as plt
import scienceplots
import numpy as np
plt.style.use('science') 

output_dir = "plots/latest_plots"

def plot_R_over_P(atm):
    '''
    Plots the pressure profile of the atmosphere as a function of radius.

    Takes input parameters: x [unit], y [unit], z [unit], ...

    All calculations done in SI units.
    '''
    plt.figure(figsize=(8,6))
    plt.plot(atm.pressures, atm.radii)
    plt.ylabel("Radius (m)")
    plt.xlabel("Pressure (Pa)")
    plt.gca().invert_xaxis() #pressure decreases with increasing radius, so we invert the x-axis to make it more intuitive
    plt.title("Pressure profile of the atmosphere")
    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.savefig(f"{output_dir}/radius_over_pressure.png", dpi=300)
    plt.close()

def plot_R_over_T(atm):
    '''
    Plots the temperature profile of the atmosphere as a function of radius.

    Takes input parameters: atmosphere object with attributes radii and temperatures.

    All calculations done in SI units.
    '''
    plt.figure(figsize=(8,6))
    plt.plot(atm.T, atm.radii)
    plt.ylabel("Radius (m)")
    plt.xlabel("Temperature (K)")
    plt.title("Temperature profile of the atmosphere")
    plt.yscale("log")
    plt.grid()
    plt.savefig(f"{output_dir}/radius_over_temperature.png", dpi=300)
    plt.close()

def plot_P_over_T(atm):
    '''
    Plots the temperature profile of the atmosphere as a function of pressure.

    Takes input parameters: atmosphere object with attributes pressures and temperatures.

    All calculations done in SI units.
    '''
    plt.figure(figsize=(8,6))
    plt.plot(atm.T, atm.pressures)
    plt.ylabel("Pressure (Pa)")
    plt.gca().invert_yaxis() #pressure decreases with increasing radius, so we invert the x-axis to make it more intuitive
    plt.xlabel("Temperature (K)")
    plt.title("Temperature profile of the atmosphere")
    plt.yscale("log")
    plt.grid()
    plt.savefig(f"{output_dir}/pressure_over_temperature.png", dpi=300)
    plt.close()

def plot_M_planet_over_M_dot(masses, M_dots, regime_break_index=None):
    '''
    Plots the mass loss rate of the atmosphere as a function of planetary mass.

    Takes input parameters: atmosphere object with attribute M_p, and mass loss rate M_dot.

    All calculations done in SI units.
    '''
    plt.figure(figsize=(8,6))
    plt.plot(M_dots, masses, label="Mass loss rate per mass")
    if regime_break_index is not None:
        plt.axhline(y=masses[regime_break_index], color='red', linestyle='--', label='Regime Break')
    plt.xlabel("Mass loss rate (kg/s)")
    plt.ylabel("Planetary mass (kg)")
    plt.title("Mass loss rate as a function of planetary mass")
    plt.grid()
    plt.legend()
    plt.savefig(f"{output_dir}/planetary_mass_over_mass_loss_rate.png", dpi=300)
    plt.close()


def plot_R_planet_over_M_dot(radii, M_dots, regime_break_index=None):
    '''
    Plots the mass loss rate of the atmosphere as a function of planetary mass.

    Takes input parameters: atmosphere object with attribute M_p, and mass loss rate M_dot.

    All calculations done in SI units.
    '''
    plt.figure(figsize=(8,6))
    plt.plot(M_dots, radii, label="Mass loss rate per radius")
    if regime_break_index is not None:
        plt.axhline(y=radii[regime_break_index], color='red', linestyle='--', label='Regime Break')
    plt.xlabel("Mass loss rate (kg/s)")
    plt.ylabel("Planetary radius (m)")
    plt.title("Mass loss rate as a function of planetary radius")
    plt.grid()
    plt.legend()
    plt.savefig(f"{output_dir}/planetary_radius_over_mass_loss_rate.png", dpi=300)
    plt.close()

def plot_xuv_flux_over_M_dot(fluxes, M_dots, regime_break_index=None, compare=False):
    '''
    Plots the mass loss rate of the atmosphere as a function of XUV flux.

    Takes input parameters: 

    All calculations done in SI units.
    '''
    plt.figure(figsize=(8,6))
    plt.plot(fluxes, M_dots, label="Modelled")
    if regime_break_index is not None:
        plt.axvline(y=fluxes[regime_break_index], color='red', linestyle='--', label='Regime Break')
    
    #comparison trend
    if compare is True:
        k = M_dots[0] / np.sqrt(fluxes[0])
        trendline = k * np.sqrt(fluxes)
        plt.plot(fluxes, trendline, label=r"Trend: $\dot{M} proportional to \sqrt{F_{XUV}}$", 
                linestyle=':', color='yellow', alpha=0.7)
    plt.ylabel("Mass loss rate (kg/s)")
    plt.xlabel("XUV flux (W m**2)")
    plt.title("Mass loss rate as a function of XUV flux received")
    plt.grid()
    plt.legend()
    plt.savefig(f"{output_dir}/xuv_flux_over_mass_loss_rate.png", dpi=300)
    plt.close()
