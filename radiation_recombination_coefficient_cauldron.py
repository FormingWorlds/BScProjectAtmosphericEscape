### This file is for calculating the radiation recombination coefficient for the different species ###

### Imports ###
import numpy as np
import pandas as pd


############################################################################

### Fetched the data manually for the species currently of interest ###
rr_coefficients_case_A = {
    # Assumes all of them to be singly ionised and T=10 000 K, and the files are all using LS coupling (no fine structure levels). Took values from the ARR(FIT) column. REFERENCE: cfout_K file from AMDPP archive, link: https://amdpp.phys.strath.ac.uk/tamoc/DATA/RR/RR_web/adf48/clist_K
    'C': 4.72 * 10**(-13), #[cm^3 s^-1]  from row Z= 6 N= 5 M= 1 W= 2        Z means nuclear charge, N means number of electrons, M means number of metastable levels, W means statistical weight of the inital state,
    'He': 4.37 * 10**(-13), #[cm^3 s^-1] from row Z= 2 N= 1 
    'O': 2.72 * 10**(-13), #[cm^3 s^-1]  from row Z= 8 N= 7 M= 1 W= 4
    'N': 3.76 * 10**(-13), #[cm^3 s^-1]  from row Z= 7 N= 6 M= 1 W= 1
}

rr_coefficients_ground_state = {
    'C':2.32 * 10**(-13), #[cm^3 s^-1] from nrb05#b_c1ls.dat, data under  PRTI= 1  TRMPRT= (2P)  SPNPRT= 2, coeff under TE index 1.00E+04
    'He': 1.56 * 10**(-13), #[cm^3 s^-1] from nrb05#h_he1ls.dat, data under PRTI= 1  TRMPRT= (2S)  SPNPRT= 2, coeff under TE index 1.00E+04
    'O': 1.31 * 10**(-13), #[cm^3 s^-1] from nrb05#n_o1ls.dat, data under  PRTI= 1  TRMPRT= (4S)  SPNPRT= 4, coeff under TE index 1.00E+04
    'N': 1.15 * 10**(-13), #[cm^3 s^-1] from nrb05#c_n1ls.dat, data under   PRTI= 1  TRMPRT= (3P)  SPNPRT= 3, coeff under TE index 1.00E+04
}


### Calculating the coefficients for the ionised species ###
rr_coefficients_case_B = {
    'H': 2.7 * 10**(-13) #[cm^3 s^-1] from Murray-Clay et. al. 2009 after eq. 7
}
for species in rr_coefficients_case_A.keys():
    rr_coefficients_case_B[species] = rr_coefficients_case_A[species] - rr_coefficients_ground_state[species] #[cm^3 s^-1] case B coefficient is case A coefficient minus the coefficient for recombinations to the ground state, where case A coefficient: radiation recombination coefficient for case A, coefficient for recombinations to the ground state: radiation recombination coefficient for recombinations to the ground state

print(f'{rr_coefficients_case_B["C"]}, {rr_coefficients_case_B["He"]}, {rr_coefficients_case_B["O"]}, {rr_coefficients_case_B["N"]}')
###########################################################################

### Must translate this to what the combined radiation recombination coefficient for the dominant species in the escaping atmosphere is, which is what we will use in the calculations. ###
#I've taken the average of the coefficients for the dominant species in the escaping atmosphere, weighted by the number of atoms of each species in the molecule, where the denominator is the total number of atoms in the molecule.
species_rr_coefficients_case_B = {
    'H': 2.7 * 10**(-13), #[cm^3 s^-1] from Murray-Clay et. al. 2009 after eq. 7
    'H2' : 2.7 * 10**(-13), #[cm^3 s^-1] from Murray-Clay et. al. 2009 after eq. 7, only valid for dissociated H2
    'H2O': (2 * rr_coefficients_case_B['H'] + rr_coefficients_case_B['O']) / 3, #[cm^3 s^-1] 
    'CO2': (rr_coefficients_case_B['C'] + 2 * rr_coefficients_case_B['O']) / 3, #[cm^3 s^-1]
    'N2': rr_coefficients_case_B['N'], #[cm^3 s^-1]
    'H2He': (2 * rr_coefficients_case_B['H'] * 0.9 + rr_coefficients_case_B['He'] * 0.1) / 3 #[cm^3 s^-1] Since it's a gas mix it should take into account ratio between gases I think.
}

print(f'{species_rr_coefficients_case_B["H2O"]}, {species_rr_coefficients_case_B["CO2"]}, {species_rr_coefficients_case_B["N2"]}, {species_rr_coefficients_case_B["H2He"]}')