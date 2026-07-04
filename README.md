# BScProjectAtmosphericEscape
This branch is part of a larger shared repository for the BSc project on atmospheric escape mechanisms (Forming Worlds Lab 2026).

In particular, this branch is dedicated to the investigation of the Radiation-Recombination-Limited (RR) hydrodynamic escape regime of exoplanets. Hydrodynamic escape is most important for smaller planets with puffy envelopes, as they are less able to hold onto their envelopes as it expands through PdV work to cool from heating (for example by stellar XUV irradiation). There are several categories of hydrodynamic escape, and RR-limited escape occurs when the atmosphere is sufficiently heated to dissociate and ionise atoms in the gas, where it proceeds to leak energy away due recombinations in the gas. As a consequence, the remaining energy to do PdV work and expand the gas is lower, thus the escape is limited by the recombinations. The aim of my work was to determine under which conditions the RR-limited regime occurs for super-Earths and sub-Neptunes. For the purpose of this thesis, "conditions" meant different planetary masses, incident XUV fluxes, and atmospheric compositions.

The code to look into the RR escape has the following structure:


# Main Escape Output files
These two files were the ones that were run to give output for the thesis. They are independent of each other.

#### plots/results_safari.ipynb
This file is the origin of every self-generated plot seen in the thesis. It sets up the species microphysics and the stellar cases to extract XUV flux from. Then it investigates and sometimes makes plot related to the following:
-Convergence test is conduction
-Parameter sweep showing escape rates and regimes for planets of different masses and compositions inspired by GJ1214b conditions (flux and temperature) 
-The impact of the ion-abundance-weighted recombination coefficient on the calculated escape rate
-Mass-radius relation for sub-Neptunes and super-Earths
-Gravitional well hydrodynamic wind criteria with respect to the implemented mass-radius relation
-Wind base pressure sensitivity on mass loss rate
-Sets up multi-species, multi-host star sweep and plots the resulting mass loss rates across the parameter space in several plots (main result of thesis)

#### road_runner.py
The main file for examining specific atmosphere. Was used on the PROTEUS data and simple atmosphere files. This file gave the output used to validate the escape model against literature.

# Escape Equations
These files build the escape calculations which are called upon in later files. Mainly plots/results_safari.ipynb

#### rr.escape.py
The calculations for the RR-limited escape are laid out in this file, and follows the formulas in Lopez (2017). The file also contains functions which allows you to input an atmosphere object and get relevant values and information about the escape. 

#### el_escape.py
This file contains the function for energy-limited escape, another hydrodynamic escape regime. Think of it as the standard one, where the escape is limited by how much energy is input to expand the atmosphere in the first place.

This file is used in rr_escape.py to do diagnostics and compare the two escape rates (imporant if you need to which which is relevant when)


# Stellar Synthesis
For my project I explored how different incident XUV and bolometric fluxes affected the escape, thus there are some files whose purpose is to create realisitc fluxes from mock stars. These files are called upon in plots/results_safari.ipynb and stellar_exploration.ipynb

#### stellar_cradle.py
This file contains the relevant functions to retrieve the XUV and bolometric the mock stars given stellar mass, age, and semi-major axis using a series of relations found in literature. 

#### stellar_exploration.ipynb
This notebook plots some of the relations using the functions of stellar_cradle.py to verify they work as they should.


# Assisting files
There files are called upon in other files that produce output relating to the escape.
#### proteus_fetch.py
This file makes my escape equations compatible with PROTEUS planetary output developed by the Forming Worlds Lab. The file is run on its own to retrieve the PROTEUS data found in atmospheres/ProteusGoogleDrive. It also contains a function that is called upon by atmospheres/atmosphere_setting.py to retrieve the volume mixing ratios of the PROTEUS atmospheres, which in turn sets the microphysics.

#### conversions.py
A small file which has some unit conversions often used and called upon for calculations, for example in plots/results_safari.ipynb and rr_escape.py.

#### bulk_param_examination.py
This file contains two functions which are highly called upon for parameter sweeps in plots/results_safari.ipynb.
It allows the creation of many simple atmosphere objects and to sweep over them and extract diagnostics about the RR and EL escape rates and which is limiting.

#### plots/plotter.py
This file contains functions which allows the user to easily plot some attributes of their atmospheres and hwo the escape rate changes with certian parameters. It is currently not called upon by any files.


# Atmosphere Construction
These files were used to calculate and create simple atmospheres which are later examined in plots/results_safari.ipynb or road_runner.py.

#### atmospheres/atmosphere_setting.py
This file was the core to be able to iterate over the parameter space and examine many possible configurations. It defines an atmosphere class such that a simple atmosphere can easily be created and where all its attributes are stored together. This file is called upon by bulk_param_examination.py, plots/results_safari.ipynb, and in the following manually configured simple atmosphere files: benchmark_planet.py, simple_H2.py, simple_H2He_mix.py, simple_H2O.py

#### benchmark_planet.py, simple_H2.py, simple_H2He_mix.py, simple_H2O.py
These files are setting up a single, simple atmosphere which can be examined. Their atmospheres must be called upon in another file for examination (except for the benchmark_planet.py, which can both be called upon and ran).

simple_H2He_mix.py also contains the atmospheres of the planets used to valide the model against Salz 2016.

#### radiation_recombination_coefficient_cauldron.py
This file calculates the species dependent ion-abundance-weighted recombination coefficient imported by atmospheres/atmosphere_setting.py. The data to calculate it was fetched manually from the folder: coefficients_calculation_data

