# BScProjectAtmosphericEscape
This repository contains the code I used in my BSc thesis titled "The Influence of Atmospheric Composition and Photodissociation on Jeans Escape from Exoplanets". It uses the PROTEUS data folder (atmospheric profiles generated with the PROTEUS framework) as input and calculates the Jeans escape rates for each given atmosphere: H₂-, H₂O-, CO₂- and N₂-dominated atmospheres for 1 and 10 Earth-mass planets and 1 and 1000 present-day Earth stellar irradiation. It identifies the exobase, performs a Bates temperature profile extension when necessary, and compares two assumptions for the upper atmosphere:

1. constant volume mixing ratios
2. complete photodissociation

The code generates all figures included in the thesis.

# Requirements

Required packages are listed in requirements.txt.

# Main Physics

The files **exobase.py**, **physics.py**, and **jeans.py** implement the core equations that are needed for the Jeans escape formalism according to Van Looveren et al. (2024), equations 1-6. 

# Running the code

The **proteus_ct_VMR.py** and **proteus_dissociation.py** are the files that should be run to generate the output CSV files of the Jeans escape rates for the two cases taken into account: constant VMR and photodissociation, respectively. They depend on **extend_ct_VMR.py** and **extend_dissociation** to create the extension in the two upper-atmosphere assumptions.
# Figure Mapping (from the Plotting Code and Validation Checks folder)

The plotting scripts read the CSV files produced by **proteus_ct_VMR.py** and **proteus_dissociation.py** to reproduce the figures in the thesis.

For the Validation section, the folder Validation Checks contains all the necessary scripts. 
1. **testcase.py** is the first check (Table 3.1)
2. **test_earth.py** is the Earth case test (Table 3.2)
3. **test_comparison.py** is the Van Looveren comparison: Figures 3.1, 3.2, 3.3

The **geometry_test.py** script tests whether the day-side area loss is exactly half of the full-side area and is just a physical check.

For the rest of the results:

1. **extension_height.py** is Figure 3.4.
2. **dominant_species_comparison.py**: Figure 3.6
3. **plotter_by_T_inf**: Figures 3.8, B.1, B.2, B.3, B.4. This code can be run for both extension cases by changing the outdir, input csv files, and the location where the figures are saved.
4. **more_plots.py**: Figures 3.5, 3.7, 3.9

The rest of the scripts in the Plotting Code were used in other stages of the thesis writing for visualisation purposes. The **pT_plots.py** creates the profile-temperature plots using the Bates extension for all cases. The **plotter.py** uses the first valid T∞ case (typically 200 K) from the summary table as a representative point, and it was used before the temperature sensitivity study was introduced.




