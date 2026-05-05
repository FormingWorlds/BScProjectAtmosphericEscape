#functions for plotting parameters of the planets and the escape stuff
import os
import matplotlib.pyplot as plt
import scienceplots
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
    plt.xscale("log")
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
    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.savefig(f"{output_dir}/pressure_over_temperature.png", dpi=300)
    plt.close()

def plot_M_planet_over_M_dot(masses, M_dots):
    '''
    Plots the mass loss rate of the atmosphere as a function of planetary mass.

    Takes input parameters: atmosphere object with attribute M_p, and mass loss rate M_dot.

    All calculations done in SI units.
    '''
    plt.figure(figsize=(8,6))
    plt.plot(M_dots, masses)
    plt.xlabel("Mass loss rate (kg/s)")
    plt.ylabel("Planetary mass (kg)")
    plt.title("Mass loss rate as a function of planetary mass")
    plt.grid()
    plt.savefig(f"{output_dir}/planetary_mass_over_mass_loss_rate.png", dpi=300)
    plt.close()

