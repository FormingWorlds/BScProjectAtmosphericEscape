# BScProjectAtmosphericEscape
repository for the BSc project on diffusion limited atmospheric escape.

This repository contains all the scripts necessary to reproduce the plots of my bachelor thesis on diffusion limited atmospheric escape. It also contains the raw PROTEUS atmosphere files that were used as input data.

see requirements.txt for the packages and versions ive used

the PROTEUS_data/ folder contains all data pertaining to the initial PROTEUS atmospheres

the analysis/ folder containts jeans escape data (and data made by proteus_output.py for convenience) which is required to run jeans vs disso plot makers


the model_validation/ folder contains all scripts pertaining to the model validation section of the thesis:
classical_limit.py is the classical limiting flux function
classical_limit_no_T_grad.py is the classical limiting flux function ignoring the temperature gradient
improved_limit.py is the improved limiting flux function
non-proteus_temp_sweep.py is the temperature sweep plotter, producing Fig.3.
non-proteus_sweep.py is the other parameter sweep plotter, producing Fig.4.
slattery_diff.py is the script calculating the slattery diffusion coefficient where needed
diffusion_coefficients.csv contains diffusion coefficient data. there is also a metadata.txt file with more information on the data.


the improved_model/ folder contains all scripts pertaining to the analysis of proteus atmospheres
the full pipeline is:
proteus_fetch.py which fetches proteus atmospheres and creates a class
proteus_extra_functions.py which contains extra needed functions
proteus_physics.py which contains the physics calculation (! this is where you change the dissociation fraction)
proteus_output.py which then saves the calculations as a .csv 

diff_vs_jeans_csv.py then creates a neat csv file for further plotting jeans vs disso and other plots. needed for plotting(!)

jeans_vs_disso_plotting.py is the script plotting figures 8-11 and 14-17, but the dissociated/non-dissociated case and species need to be changed manually (lines 60,63,70,78,85,93,100,108,115,134,137,140,143. either string of species name, or lim/jeans)

proteus_output_profile_plotting.py is the script plotting all 16 profile analyses at once, as the one in Fig. 2
proteus_physics_for_profiles.py is the script needed for the above plotter to work (just a cut of proteus_physics.py)

proteus_physics_for_D_vs_T.py is the script needed for the plotter of D vs T to work (contained in plot_boogaloo.ipynb)

plot_boogaloo.ipynb is a notebook containing all the other plots not mentioned that were used (fig 5,6,7,12,13). first line is a comment on which plot it is (which fig number). some require specific files which are made by one of the above menitioned scripts. if run in order, will work.

slattery_diff.py is the script calculating the slattery diffusion coefficient where needed
