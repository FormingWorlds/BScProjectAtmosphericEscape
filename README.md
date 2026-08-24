# BScProjectAtmosphericEscape
This branch is part of the shared repository BScProjectAtmosphericEscape (2026), part of the Forming Worlds Lab.

This repository contains the code used in the Bachelor Thesis named 'Impact-driven Atmospheric Escape for Various Impactor Sizes and Size Distributions' by Renske Beuker. It focuses on modeling the atmospheric escape caused by impacts. This is based on the paper by Schlichting et al. (2015): Atmospheric mass loss during planet formation: The importance of planetesimal impacts. The code here builds upon this paper, researching differing values for mass loss rate.

# Requirements
The required versions of packages used are given in requirements.txt.

# PROTEUS_data/
The folder <b>PROTEUS_data/</b> contains the model atmospheres simulated by the PROTEUS framework, specifically AGNI. The results are based on these atmospheres. The model atmospheres are made for planet masses of 1 and 10 M_Earth, 1 and 1000 F_Earth, and H2, H2O, N2, and CO2.

The data is accessed using the file <b>open_PROTEUS_csv.py</b>. This file contains functions using pandas that read the files and add the needed information to a dictionary.

# The Main Physics
<b>properties.py</b> contains the functions to calculate scale height of the atmosphere, as well as the total atmosphere mass. The total atmosphere mass is used for global mass loss (<b>giant_impacts_fractional_loss.py</b>), by multiplying the fractional loss with the atmosphere mass. The scale height is used for planetesimal mass loss, as it determines radius thresholds for the impactor as well as cap mass (<b>planetesimal_and_MLR.py</b>).

<b>giant_impacts_fractional_loss.py</b> contains the physics for global mass loss caused by giant impacts. An adiabatic atmosphere is used to model this. The outputs of this code are found in <b>Outputs/Global_mass_loss_data.csv</b>. It contains the mass lost in kg for the radius of the impactor (m), as well as for the ratio m_impactor/M_planet. It also contains atmospheric mass.


<b>planetesimal_and_MLR.py</b> contains a lot of physics. It contains the radius thresholds: r_min, r_cap, r_gi. These functions are used to determine the needed size of the impactor. r_min is the minimum size an impactor needs to be to eject mass, r_gi is the size where global mass loss can happen, and r_cap is the radius at which the total cap mass is lost. The cap mass function is also given in this file, as well as the ejected mass for a single planetesimal. The ejected mass by one planetesimal impactor is given as output in <b>Outputs/Planetesimal_mass_loss_data.csv</b>, together with radius thresholds, cap mass, and atmospheric mass.

The mass loss rate (MLR) for a distribution of impactors is also calculated in <b>planetesimal_and_MLR.py</b>. MLR is compared to differential power law index q, and to total impactor mass per unit time M_pl. The outputs for these are given in csv-files, under Outputs/Mass_loss_rate_vs_q.csv and <b>Outputs/Mass_loss_rate_vs_M_pl.csv</b>. These files also contain cap mass and radius thresholds.

<b>M_imp_tot_per_rad.py</b> gives a function to recreate figure 14 by Schlichting et al., showing how much impactor mass is needed for different impactor radii to eject the full atmosphere. The results for the model atmospheres are given in <b>Plot_code/fig14_models.py</b>.

# Running the Code
The files mentioned above contain the main physics and give the outputs. Running said files would give the output files. As they are already given, this is no longer necessary to run. However, they can be run by changing the values for inputs, if one wants to research what effect a change in variables can have. The files under <b>Plot_code/</b> need to be run to produce the plots. These files use the date from <b>Outputs/</b> and produce figures. The figures are presented in <b>figures.ipynb</b>, so that this is not necessary.

# Outputs/
Like mentioned above, outputs from the code are given in the folder <b>Outputs/</b>.
<b>Outputs/Atmospheric_profiles</b> includes the given PROTEUS profiles, along with planet and atmosphere mass, as well as effective scale height h_eff.

# Plotting figures
The code for plotting figures is given under <b>Plot_code/</b>. The figures themselves can be seen in the notebook <b>figures.ipynb</b>, so there is no need to run the <b>Plot_code/</b>.

# Validation/
The folder <b>Validation/</b> contains <b>frac_global_adia_vs_iso.py</b>, a code to plot the difference between isothermal and adiabatic atmospheres. It also contains <b>radius_cutoffs_validation.py</b>. Running this will show the radius thresholds for Earth conditions. They agree with Schlichting et al. (2015). Finally, it has <b>Earth_total_impactor_mass_validation.py</b>, to recreate figure 14 by Schlichting et al. and see how much mass is needed to eject total atmospheres for different radii.
