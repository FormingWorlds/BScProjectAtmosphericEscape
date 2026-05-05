"""This file makes a class that can hold atmosphere parameters such that they don't have to be input as separate arguments to the rr_escape_rate function. 
This is useful for testing the function with different atmospheric compositions, and for making it easier to use the function in a more general context."""
import numpy as np

class Atmosphere:
    def __init__(self, T_wind, mu_wind, M_p, F_ins, nu_0, mu_plus_wind, R_p, pressures, temperatures, heights, vmrs=None, P_base=10**(-4), microphysics_cap=None):

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
        self.microphysics_cap = microphysics_cap if microphysics_cap is not None else 0.01

    def calc_xuv_from_insolation(self, F_ins):
        '''
        Calculates the XUV flux from the total insolation.

        Takes input parameters: x [unit], y [unit], z [unit], ...

        All calculations done in SI units.
        '''
        xuv_fraction = 10**(-6) #REFERENCE Murray-Clay et. al. 2009, XUV flux is typically 10^-6 times the total insolation for a sunlike star
        F_xuv = F_ins * xuv_fraction #take XUV fraction of the total insolation. Assumes sunlike star. From Murray-Clay et. al. 2009

        return F_xuv

    def read_off_wind_base_parameters(self):
        '''
        Reads off the necessary parameters for the radiation-recombination-limited escape rate calculation at the base of the escaping atmosphere.

        Takes input parameters: x [unit], y [unit], z [unit], ...

        All calculations done in SI units.
        '''
        #finds the absolute difference between the pressure profile and the target pressure at the base of the escaping atmosphere
        difference_array = np.absolute(self.pressures - self.P_base) #[Pa] array of the absolute difference between the pressure profile and the pressure at the base of the escaping atmosphere, where pressures: pressure profile of the atmosphere based on the barometric formula, P_base: pressure at the base of the escaping atmosphere

        #finds the index of minimum element from the array
        index = difference_array.argmin()
        self.R_base = self.radii[index] #[m] radius of the base of the escaping atmosphere, where radii: array of radii from the planetary radius to 10 times the planetary radius, index: index of minimum element from the array of the absolute difference between the pressure profile and the pressure at the base of the escaping atmosphere
        self.T_base = self.T[index] #[K] temperature at the base of the escaping atmosphere, where self.T: array of temperatures as a function of radius, index: index of the radius of the base of the escaping atmosphere
        
        #find vmrs at base
        self.vmrs_base = {species: vmr_array[index] for species, vmr_array in self.vmrs.items()} #vmrs at the base of the escaping atmosphere, where self.vmrs: dictonary of arrays of volume mixing ratios as a function of radius, index: index of the radius of the base of the escaping atmosphere
        #find dominant species at base
        self.dominant_species = max(self.vmrs_base, key=self.vmrs_base.get)
    

    #def determine_wind_microphysics(self):
        #######if self.microphysics_cap is None:
            #relevant_species = {species: vmr for species, vmr in self.vmrs_base.items() if vmr > self.microphysics_cap} #finds the relevant species at the base of the escaping atmosphere, where self.vmrs_base: dictonary of volume mixing ratios at the base of the escaping atmosphere, cap: threshold for determining whether a species is relevant or not based on its volume mixing ratio
        


        

        


    