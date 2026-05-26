### This file is for calculating the radiation recombination coefficient for the different species ###

import numpy as np

def calculate_rr_coefficient_case_A(T, T_0, T_1, A, B, C=None, T_2=None):
    '''
    Calculates the radiation recombination coefficient for case A for your species and environment of interest.
    Case A means that all recombinations are included, in particular the ground state is included, so that recombinations to the ground state are not followed by immediate re-ionisation. 

    Takes input parameters: T is the temperature [K], T_0 [K], T_1 [K], A [cm^3 s^-1], B [dimensionless], optional parameters C [dimensionless], T_2 [K] for low charge atoms.

    Outputs: radiation recombination coefficient [cm^3 s^-1]

    Designed to take input parameters from data from AMDPP (Atomic and Molecular Diagnostics Processes in Plasma).
    '''

    # Checks if the optional parameters C and T_2 are provided, meaning we are working with a low charge atom, and if so, calculates the effective B parameter (B_eff) using the provided formula. If not, it sets B_eff to be equal to B.
    if C is not None and T_2 is not None:
        B_eff = B + C * np.exp(-T_2 / T)
    else:
        B_eff = B

    # Calculates the radiation recombination coefficient for case A using the provided formula, which incorporates the temperature dependence and the effective B parameter.
    # REFERENCE AMDPP Fits for Total Radiative Recombination Coefficients, link: https://amdpp.phys.strath.ac.uk/tamoc/DATA/RR/RR_web/adf48/
    alpha_rr_case_A = A * ( np.sqrt(T / T_0) * (1 + np.sqrt(T / T_0))**(1-B_eff) * (1 + np.sqrt(T / T_1))**(1+B_eff) )**(-1)

    return alpha_rr_case_A