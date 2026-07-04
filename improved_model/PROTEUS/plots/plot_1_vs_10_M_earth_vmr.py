import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

atm_archetype = ["CO2", "H2", "H2O", "N2"]  
instellation = ["1_F_earth", "1000_F_earth"]
mass = ["1_M_earth", "10_M_earth"]

mass_labels = {
    "1_M_earth": r"$M_p = 1\,M_{\oplus}$",
    "10_M_earth": r"$M_p = 10\,M_{\oplus}$"
}

# Define consistent colorblind-friendly colors for each species
# Using Paul Tol's colorblind-safe palette
species_colors = {
    # Atmosphere archetypes (distinctive colors)
    "CO2": "#D55E00",   # Vermillion (orange-red)
    "H2": "#0072B2",    # Blue
    "H2O": "#56B4E9",   # Sky blue
    "N2": "#CC79A7",    # Reddish purple
    
    # Other major species
    "H": "#E69F00",     # Orange
    "O": "#009E73",     # Bluish green
    "C": "#F0E442",     # Yellow
    "CH4": "#882255",   # Purple
    "CO": "#AA4499",    # Mauve
    "O2": "#44AA99",    # Teal
    "N": "#999933",     # Olive
    "NH3": "#117733",   # Green
    "H2S": "#88CCEE",   # Cyan
    "S": "#DDCC77",     # Sand
    "S2": "#661100",    # Dark red
    "SO2": "#332288"    # Indigo
}

# Plot species vmr vs pressure, if does not find file print warning and skip
species_vmr = ["H2","H2O","H","O","C","CH4","CO","CO2","O2","N","N2","NH3","H2S","S","S2","SO2"]


# Plot VMR vs pressure for each atm_archetype
for sp in atm_archetype:
    # Create figure with 2 panels side by side
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    plotted_species = set()  # Track which species have been plotted
    # Loop over mass (one panel per mass)
    for idx, mas in enumerate(mass):
        ax = axes[idx]
        for inst in instellation:
            input_dir = f"../{sp}_atmospheres/{inst}/"
            csv_path = os.path.join(input_dir, f'{sp}_atmosphere_{mas}_{inst}.csv')
            if not os.path.exists(csv_path):
                print(f"File not found: {csv_path}")
                continue
            df = pd.read_csv(csv_path, sep="\t")
            
            # Plot all species for this instellation
            for s_vmr in species_vmr:
                col = s_vmr + " [VMR]"
                if col in df.columns:
                    label = s_vmr if s_vmr not in plotted_species else None

                    if inst == "1_F_earth":
                        if s_vmr == f"{sp}":
                            ax.plot(df[col], df["Pressure [Pa]"]*1e-5, label=label, linestyle='-', linewidth=3, color=species_colors.get(s_vmr, None))
                        else:
                            ax.plot(df[col], df["Pressure [Pa]"]*1e-5, label=label, linestyle='-', linewidth=1.5, color=species_colors.get(s_vmr, None))
                    else:
                        if s_vmr == f"{sp}":
                            ax.plot(df[col], df["Pressure [Pa]"]*1e-5, label=label, linestyle='--', linewidth=3, color=species_colors.get(s_vmr, None))
                        else:
                            ax.plot(df[col], df["Pressure [Pa]"]*1e-5, label=label, linestyle='--', linewidth=1.5, color=species_colors.get(s_vmr, None))
                    plotted_species.add(s_vmr)

        # Create custom legend for instellation (line styles)
        instellation_legend_elements = [
            Line2D([0], [0], color='black', linestyle='-', lw=2, label=r'1 $F_{\oplus}$'),
            Line2D([0], [0], color='black', linestyle='--', lw=2, label=r'1000 $F_{\oplus}$')
        ]   

        species_handles = [
            Line2D([0], [0], color=species_colors[s], lw=2, label=s)
            for s in species_vmr
            if s in species_colors
        ]

        # Add both legends
        first_legend = ax.legend(handles=instellation_legend_elements, loc='lower center', title='Instellation')
        ax.add_artist(first_legend)  # Keep first legend when adding second

        species_legend = ax.legend(
            handles=species_handles,
            loc='lower left',
            title="Species",
            fontsize=10
        )
        ax.add_artist(species_legend)

        ax.set_yscale("log")
        ax.invert_yaxis()
        ax.set_xscale("log")
        ax.tick_params(axis='both', which='major', labelsize=14)
        ax.set_xlim(1e-10, 2)
        ax.set_ylim(3e4, 0.5e-10)
        ax.set_xlabel("Volume Mixing Ratio", fontsize=14)
        ax.set_ylabel("Pressure [bar]", fontsize=14)
        ax.get_legend()
        ax.set_title(mass_labels[mas], fontsize=14)

    plt.suptitle(f'{sp}-dominated atmospheres', fontsize=14)
    plt.savefig(f"VMR/VMR_vs_Pressure_{sp}_1_vs_10_M_earth.png", dpi= 300, bbox_inches='tight')
    plt.close()