import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from netCDF4 import Dataset

species = ["CO2", "H2", "H2O", "N2"]
instellation = ["1_F_earth", "1000_F_earth"]
mass = ["1_M_earth", "10_M_earth"]

species_colors = {
    "CO2": "#D55E00",
    "H2": "#0072B2",
    "H2O": "#009E73",
    "N2": "#CC79A7"
}

def extract_atm_columns(nc_path):
    with Dataset(nc_path, 'r') as ds:
        pressure = ds.variables['pl'][:]
        temperature = ds.variables['tmpl'][:]
        kzz = ds.variables['Kzz'][:]
        
        if pressure.ndim == 2:
            pressure = pressure[-1, :]
            temperature = temperature[-1, :]
            kzz = kzz[-1, :]
        
        return pressure, temperature, kzz

##########################################
# Plot temperature vs pressure
##########################################

# Create figure with 2 panels side by side
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Loop over mass (one panel per mass)
for idx, mas in enumerate(mass):
    ax = axes[idx]
    
    for inst in instellation:
        for sp in species:
            input_dir = f"../{sp}_atmospheres/{inst}/"
            csv_path = os.path.join(input_dir, f'{sp}_atmosphere_{mas}_{inst}.csv')
            
            if not os.path.exists(csv_path):
                print(f"Missing: {csv_path}")
                continue
            
            df = pd.read_csv(csv_path, sep='\t')
            pressure = df['Pressure [Pa]'] * 1e-5  # Pa → bar
            temperature = df['Temperature [K]']

            # Line style based on INSTELLATION
            linestyle = '-' if inst == "1_F_earth" else '--'
            
            ax.plot(
                temperature,
                pressure,
                linestyle=linestyle,
                linewidth=2,
                color=species_colors[sp]
            )
    
    # --- LEGENDS for each panel ---
    species_legend = [
        Line2D([0], [0], color=species_colors[sp], lw=2, label=sp)
        for sp in species
    ]
    
    instellation_legend = [
        Line2D([0], [0], color='black', linestyle='-', lw=2, label=r'1 $F_{\oplus}$'),
        Line2D([0], [0], color='black', linestyle='--', lw=2, label=r'1000 $F_{\oplus}$')
    ]
    
    # Add both legends
    first_legend = ax.legend(handles=species_legend, loc='upper right', title='Species')
    ax.add_artist(first_legend)
    ax.legend(handles=instellation_legend, loc='upper center', title='Instellation')
    
    # --- AXES ---
    ax.set_yscale('log')
    ax.invert_yaxis()
    # ax.set_ylim(1e-9, 0.5e-10)
    # ax.set_xlim(0, 1000)
    ax.tick_params(labelsize=14)
    ax.set_xlabel('Temperature [K]', fontsize=14)
    ax.set_ylabel('Pressure [bar]', fontsize=14)
    
    # Set title based on mass
    if mas == "1_M_earth":
        ax.set_title(r'$M_{planet} = 1\,M_{\oplus}$', fontsize=14)
    else:
        ax.set_title(r'$M_{planet} = 10\,M_{\oplus}$', fontsize=14)

plt.tight_layout()
plt.savefig('profiles/Temperature_vs_Pressure_1_vs_10_M_earth.png', dpi=300)
plt.close()

##########################################
# Plot temperature vs height
##########################################

# Create figure with 2 panels side by side
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
# Loop over mass (one panel per mass)
for idx, mas in enumerate(mass):
    ax = axes[idx]
    
    for inst in instellation:
        for sp in species:
            input_dir = f"../{sp}_atmospheres/{inst}/"
            csv_path = os.path.join(input_dir, f'{sp}_atmosphere_{mas}_{inst}.csv')
            
            if not os.path.exists(csv_path):
                print(f"Missing: {csv_path}")
                continue
            
            df = pd.read_csv(csv_path, sep='\t')
            height = df['Height [m]']
            temperature = df['Temperature [K]']

            # Line style based on INSTELLATION
            linestyle = '-' if inst == "1_F_earth" else '--'
            
            ax.plot(
                temperature,
                height,
                linestyle=linestyle,
                linewidth=2,
                color=species_colors[sp]
            )
    
    # --- LEGENDS for each panel ---
    species_legend = [
        Line2D([0], [0], color=species_colors[sp], lw=2, label=sp)
        for sp in species
    ]
    
    instellation_legend = [
        Line2D([0], [0], color='black', linestyle='-', lw=2, label=r'1 $F_{\oplus}$'),
        Line2D([0], [0], color='black', linestyle='--', lw=2, label=r'1000 $F_{\oplus}$')
    ]
    
    # Add both legends
    first_legend = ax.legend(handles=species_legend, loc='upper right', title='Species')
    ax.add_artist(first_legend)
    ax.legend(handles=instellation_legend, loc='lower left', title='Instellation')
    
    # --- AXES ---
    ax.set_yscale('symlog')
    #ax.invert_yaxis()
    # ax.set_ylim(1e-9, 0.5e-10)
    ax.tick_params(labelsize=14)
    ax.set_xlabel('Temperature [K]', fontsize=14)
    ax.set_ylabel('Height [m]', fontsize=14)
    
    # Set title based on mass
    if mas == "1_M_earth":
        ax.set_title(r'$M_{planet} = 1\,M_{\oplus}$', fontsize=14)
    else:
        ax.set_title(r'$M_{planet} = 10\,M_{\oplus}$', fontsize=14)

plt.tight_layout()
plt.savefig('profiles/Temperature_vs_Height_1_vs_10_M_earth.png', dpi=300)
plt.close()

##########################################
# Plot Kzz vs pressure
##########################################
# Create figure with 2 panels side by side
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
# Loop over mass (one panel per mass)
for idx, mas in enumerate(mass):
    ax = axes[idx]
    
    for inst in instellation:
        for sp in species:
            input_dir = f"../{sp}_atmospheres/{inst}/"
            csv_path = os.path.join(input_dir, f'{sp}_atmosphere_{mas}_{inst}.csv')
            
            if not os.path.exists(csv_path):
                print(f"Missing: {csv_path}")
                continue
            
            df = pd.read_csv(csv_path, sep='\t')
            pressure = df['Pressure [Pa]'] * 1e-5  # Pa → bar
            Kzz = df['Kzz [cm2/s]']

            # Line style based on INSTELLATION
            linestyle = '-' if inst == "1_F_earth" else '--'
            
            ax.plot(
                Kzz,
                pressure,
                linestyle=linestyle,
                linewidth=2,
                color=species_colors[sp]
            )
    
    # --- LEGENDS for each panel ---
    species_legend = [
        Line2D([0], [0], color=species_colors[sp], lw=2, label=sp)
        for sp in species
    ]
    
    instellation_legend = [
        Line2D([0], [0], color='black', linestyle='-', lw=2, label=r'1 $F_{\oplus}$'),
        Line2D([0], [0], color='black', linestyle='--', lw=2, label=r'1000 $F_{\oplus}$')
    ]
    
    # Add both legends
    first_legend = ax.legend(handles=species_legend, loc='lower right', title='Species')
    ax.add_artist(first_legend)
    ax.legend(handles=instellation_legend, loc='upper left', title='Instellation')
    
    # --- AXES ---
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.invert_yaxis()
    #ax.set_ylim(3e4, 0.5e-10)
    ax.tick_params(labelsize=14)
    ax.set_xlabel(r'Kzz [cm$^2$/s]', fontsize=14)
    ax.set_ylabel('Pressure [bar]', fontsize=14)
    
    # Set title based on mass
    if mas == "1_M_earth":
        ax.set_title(r'$M_{planet} = 1\,M_{\oplus}$', fontsize=14)
    else:
        ax.set_title(r'$M_{planet} = 10\,M_{\oplus}$', fontsize=14)

plt.tight_layout()
plt.savefig('profiles/Kzz_vs_Pressure_1_vs_10_M_earth.png', dpi=300)
plt.close()

##########################################
# Plot density vs pressure
##########################################

# Create figure with 2 panels side by side
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Loop over mass (one panel per mass)
for idx, mas in enumerate(mass):
    ax = axes[idx]
    
    for inst in instellation:
        for sp in species:
            input_dir = f"../{sp}_atmospheres/{inst}/"
            csv_path = os.path.join(input_dir, f'{sp}_atmosphere_{mas}_{inst}.csv')
            
            if not os.path.exists(csv_path):
                print(f"Missing: {csv_path}")
                continue
            
            df = pd.read_csv(csv_path, sep='\t')
            pressure = df['Pressure [Pa]'] * 1e-5  # Pa → bar
            density = df['Density [kg/m3]']

            # Line style based on INSTELLATION
            linestyle = '-' if inst == "1_F_earth" else '--'
            
            ax.plot(
                density,
                pressure,
                linestyle=linestyle,
                linewidth=2,
                color=species_colors[sp]
            )
    
    # --- LEGENDS for each panel ---
    species_legend = [
        Line2D([0], [0], color=species_colors[sp], lw=2, label=sp)
        for sp in species
    ]
    
    instellation_legend = [
        Line2D([0], [0], color='black', linestyle='-', lw=2, label=r'1 $F_{\oplus}$'),
        Line2D([0], [0], color='black', linestyle='--', lw=2, label=r'1000 $F_{\oplus}$')
    ]
    
    # Add both legends
    first_legend = ax.legend(handles=species_legend, loc='upper right', title='Species')
    ax.add_artist(first_legend)
    ax.legend(handles=instellation_legend, loc='upper center', title='Instellation')
    
    # --- AXES ---
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.invert_yaxis()
    #ax.set_ylim(3e4, 0.5e-10)
    ax.tick_params(labelsize=14)
    ax.set_xlabel(r'Density [kg/m$^3$]', fontsize=14)
    ax.set_ylabel('Pressure [bar]', fontsize=14)
    
    # Set title based on mass
    if mas == "1_M_earth":
        ax.set_title(r'$M_{planet} = 1\,M_{\oplus}$', fontsize=14)
    else:
        ax.set_title(r'$M_{planet} = 10\,M_{\oplus}$', fontsize=14)

plt.tight_layout()
plt.savefig('profiles/Density_vs_Pressure_1_vs_10_M_earth.png', dpi=300)
plt.close()

##########################################
# Plot MMW vs pressure
##########################################

# Create figure with 2 panels side by side
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Loop over mass (one panel per mass)
for idx, mas in enumerate(mass):
    ax = axes[idx]
    
    for inst in instellation:
        for sp in species:
            input_dir = f"../{sp}_atmospheres/{inst}/"
            csv_path = os.path.join(input_dir, f'{sp}_atmosphere_{mas}_{inst}.csv')
            
            if not os.path.exists(csv_path):
                print(f"Missing: {csv_path}")
                continue
            
            df = pd.read_csv(csv_path, sep='\t')
            pressure = df['Pressure [Pa]'] * 1e-5  # Pa → bar
            density = df['MMW [g/mol]']

            # Line style based on INSTELLATION
            linestyle = '-' if inst == "1_F_earth" else '--'
            
            ax.plot(
                density,
                pressure,
                linestyle=linestyle,
                linewidth=2,
                color=species_colors[sp]
            )
    
    # --- LEGENDS for each panel ---
    species_legend = [
        Line2D([0], [0], color=species_colors[sp], lw=2, label=sp)
        for sp in species
    ]
    
    instellation_legend = [
        Line2D([0], [0], color='black', linestyle='-', lw=2, label=r'1 $F_{\oplus}$'),
        Line2D([0], [0], color='black', linestyle='--', lw=2, label=r'1000 $F_{\oplus}$')
    ]
    
    # Add both legends
    first_legend = ax.legend(handles=species_legend, loc='lower left', title='Species')
    ax.add_artist(first_legend)
    ax.legend(handles=instellation_legend, loc='lower center', title='Instellation')
    
    # --- AXES ---
    ax.set_yscale('log')
    #ax.set_xscale('log')
    ax.invert_yaxis()
    #ax.set_ylim(3e4, 0.5e-10)
    ax.tick_params(labelsize=14)
    ax.set_xlabel('Mean molecular weight [g/mol]', fontsize=14)
    ax.set_ylabel('Pressure [bar]', fontsize=14)
    
    # Set title based on mass
    if mas == "1_M_earth":
        ax.set_title(r'$M_{planet} = 1\,M_{\oplus}$', fontsize=14)
    else:
        ax.set_title(r'$M_{planet} = 10\,M_{\oplus}$', fontsize=14)

plt.tight_layout()
plt.savefig('profiles/MMW_vs_Pressure_1_vs_10_M_earth.png', dpi=300)
plt.close()

##########################################
# Plot WL(nm) vs Flux(ergs/cm**2/s/nm) from star.csv file
##########################################
input_dir = f"../"
csv_path = os.path.join(input_dir, 'star.dat')
if not os.path.exists(csv_path):
    print(f"File not found: {csv_path}")
else:
    star_df = pd.read_csv(csv_path, sep=r'\s+', comment='#', names=['WL [nm]', 'Flux [ergs/cm^2/s/nm]'])
    plt.plot(star_df['WL [nm]'], star_df['Flux [ergs/cm^2/s/nm]'], 
                color='black', 
                linestyle='-',
                label='Sun',
                linewidth=2)

plt.xscale('log')
plt.yscale('log')
# plt.xlim(1, 1e4)
# plt.ylim(1e-3, 1e4)
plt.legend(fontsize=12, loc='best')
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Wavelength [nm]', fontsize=14)
plt.ylabel(r'Flux [ergs/cm$^2$/s/nm]', fontsize=14)
plt.tight_layout()
plt.savefig('profiles/Stellar_Flux.png', dpi=300)