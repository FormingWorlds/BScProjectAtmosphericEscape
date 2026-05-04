#functions for plotting parameters of the planets and the escape stuff
import matplotlib.pyplot as plt
import scienceplots
plt.style.use('science') 

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
    plt.show()


