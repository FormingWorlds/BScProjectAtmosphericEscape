# BScProjectAtmosphericEscape
This branch is part of the shared repository BScProjectAtmosphericEscape (2026), part of the Forming Worlds Lab.

This repository contains the code used in the Bachelor Thesis named 'Impact-driven Atmospheric Escape for Various Impactor Sizes and Size Distributions' by Renske Beuker. It focuses on modeling the atmospheric escape caused by impacts. This is based on the paper by Schlichting et al. (2015): Atmospheric mass loss during planet formation: The importance of planetesimal impacts. The code here builds upon this paper, researching differing values for mass loss rate.

# Requirements
The required versions of packages used are given in requirements.txt.

# PROTEUS_data/
The folder PROTEUS_data/ contains the model atmospheres simulated by the PROTEUS framework, specifically AGNI. The results are based on these atmospheres. The model atmospheres are made for planet masses of 1 and 10 M_Earth, 1 and 1000 F_Earth, and H2, H2O, N2, and CO2.

The data is accessed using the file open_PROTEUS_csv.py. This file contains functions using pandas that read the files and add the needed information to a dictionary.

# The Main Physics
<b>properties.py</b> contains the functions to calculate scale height of the atmosphere, as well as the total atmosphere mass. The total atmosphere mass is used for global mass loss (giant_impacts_fractional_loss.py), by multiplying the fractional loss with the atmosphere mass. The scale height is used for planetesimal mass loss, as it determines radius thresholds for the impactor as well as cap mass (planetesimal_and_MLR.py).

giant_impacts_fractional_loss.py contains the physics for global mass loss caused by giant impacts. An adiabatic atmosphere is used to model this. The outputs of this code are found in Outputs/Global_mass_loss_data.csv. It contains the mass lost in kg for the radius of the impactor (m), as well as for the ratio m_impactor/M_planet. It also contains atmospheric mass.


planetesimal_and_MLR.py contains a lot of physics. It contains the radius thresholds: r_min, r_cap, r_gi. These functions are used to determine the needed size of the impactor. r_min is the minimum size an impactor needs to be to eject mass, r_gi is the size where global mass loss can happen, and r_cap is the radius at which the total cap mass is lost. The cap mass function is also given in this file, as well as the ejected mass for a single planetesimal. The ejected mass by one planetesimal impactor is given as output in Outputs/Planetesimal_mass_loss_data.csv, together with effective scale height, radius thresholds, cap mass, and atmospheric mass.

The mass loss rate (MLR) for a distribution of impactors is also calculated in planetesimal_and_MLR.py. MLR is compared to differential power law index q, and to total impactor mass per unit time M_pl. The outputs for these are given in csv-files, under Outputs/Mass_loss_rate_vs_q.csv and Outputs/Mass_loss_rate_vs_M_pl.csv. These files also contain cap mass and radius thresholds.

M_imp_tot_per_rad.py gives a function to recreate figure 14 by Schlichting et al., showing how much impactor mass is needed for different impactor radii to eject the full atmosphere. The results for the model atmospheres are given in Plot_code/fig14_models.py.

# Running the Code
The files mentioned above contain the main physics and give the outputs. Running said files would give the output files. As they are already given, this is no longer necessary to run. However, they can be run by changing the values for inputs, if one wants to research what effect a change in variables can have. The files under Plot_code/ need to be run to produce the plots. These files use the date from Outputs/ and produce figures. The figures are presented in figures.ipynb, so that this is not necessary.

# Outputs/
Like mentioned above, outputs from the code are given in the folder Outputs/.

# Plotting figures
The code for plotting figures is given under Plot_code/. The figures themselves can be seen in the notebook figures.ipynb, so there is no need to run the Plot_code/.

# Validation/
The folder Validation/ contains frac_global_adia_vs_iso.py, a code to plot the difference between isothermal and adiabatic atmospheres. It also contains radius_cutoffs_validation.py. Running this will show the radius thresholds for Earth conditions. They agree with Schlichting et al. (2015). Finally, it has Earth_total_impactor_mass_validation.py, to recreate figure 14 by Schlichting et al. and see how much mass is needed to eject total atmospheres for different radii.
