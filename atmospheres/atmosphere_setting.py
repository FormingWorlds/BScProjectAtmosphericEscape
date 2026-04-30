"""This file makes a class that can hold atmosphere parameters such that they don't have to be input as separate arguments to the rr_escape_rate function. 
This is useful for testing the function with different atmospheric compositions, and for making it easier to use the function in a more general context."""


class Atmosphere:
    def __init__(self, T_wind, mu_wind, M_p, F_xuv, nu_0, mu_plus_wind, radii, pressures, temperatures, heights):

        #from planet bulk proerties
        self.M_p = M_p
        self.R_p = R_p
        self.F_xuv = F_xuv

        #composition dependent properties
        self.T_wind = T_wind
        self.nu_0 = nu_0
        self.mu_wind = mu_wind 
        self.mu_plus_wind = mu_plus_wind

        #imported atmospheric profiles of different properties as a function of radius
        self.T = temperatures
        self.radii = M_p + heights
        self.pressures = pressures