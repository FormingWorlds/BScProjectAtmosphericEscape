#functions for plotting parameters of the planets and the escape stuff
import os
import matplotlib.pyplot as plt
import scienceplots
plt.style.use('science') 

output_dir = "plots/latest_plots"

def plot_P_over_R(atm):
    '''
    Plots the pressure profile of the atmosphere as a function of radius.

    Takes input parameters: x [unit], y [unit], z [unit], ...

    All calculations done in SI units.
    '''
    plt.figure(figsize=(8,6))
    plt.plot(atm.radii, atm.pressures)
    plt.xlabel("Radius (m)")
    plt.ylabel("Pressure (Pa)")
    plt.title("Pressure profile of the atmosphere")
    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.savefig(f"{output_dir}/pressure_over_radius.png", dpi=300)
    plt.close()

def plot_T_over_R(atm):
    '''
    Plots the temperature profile of the atmosphere as a function of radius.

    Takes input parameters: atmosphere object with attributes radii and temperatures.

    All calculations done in SI units.
    '''
    plt.figure(figsize=(8,6))
    plt.plot(atm.radii, atm.T)
    plt.xlabel("Radius (m)")
    plt.ylabel("Temperature (K)")
    plt.title("Temperature profile of the atmosphere")
    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.savefig(f"{output_dir}/temperature_over_radius.png", dpi=300)
    plt.close()

def plot_T_over_P(atm):
    '''
    Plots the temperature profile of the atmosphere as a function of pressure.

    Takes input parameters: atmosphere object with attributes pressures and temperatures.

    All calculations done in SI units.
    '''
    plt.figure(figsize=(8,6))
    plt.plot(atm.pressures, atm.T)
    plt.xlabel("Pressure (Pa)")
    plt.ylabel("Temperature (K)")
    plt.title("Temperature profile of the atmosphere")
    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.savefig(f"{output_dir}/temperature_over_pressure.png", dpi=300)
    plt.close()

