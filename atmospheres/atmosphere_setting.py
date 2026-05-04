"""This file makes a class that can hold atmosphere parameters such that they don't have to be input as separate arguments to the rr_escape_rate function. 
This is useful for testing the function with different atmospheric compositions, and for making it easier to use the function in a more general context."""


class Atmosphere:
    def __init__(self, T_wind, mu_wind, M_p, F_ins, nu_0, mu_plus_wind, R_p, pressures, temperatures, heights, vmrs=None):

        #from planet bulk proerties
        self.M_p = M_p
        self.R_p = R_p
        self.F_xuv = self.calc_xuv_from_insolation(F_ins)

        #composition dependent properties, might need to add some function which calculates these based on what is found in the vmrs
        self.T_wind = T_wind
        self.nu_0 = nu_0
        self.mu_wind = mu_wind 
        self.mu_plus_wind = mu_plus_wind

        #imported atmospheric profiles of different properties as a function of radius
        self.T = temperatures
        self.radii = R_p + heights
        self.pressures = pressures
        self.vmrs = vmrs if vmrs is not None else {} #must be a dictonary even if no vmrs are given so that functions later don't throw a fit.

    def calc_xuv_from_insolation(self, F_ins):
        '''
        Calculates the XUV flux from the total insolation.

        Takes input parameters: x [unit], y [unit], z [unit], ...

        All calculations done in SI units.
        '''
        xuv_fraction = 10**(-6) #REFERENCE Murray-Clay et. al. 2009, XUV flux is typically 10^-6 times the total insolation for a sunlike star
        F_xuv = F_ins * xuv_fraction #take XUV fraction of the total insolation. Assumes sunlike star. From Murray-Clay et. al. 2009

        return F_xuv