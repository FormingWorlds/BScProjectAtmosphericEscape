"""This file makes a class that can hold atmosphere parameters such that they don't have to be input as separate arguments to the rr_escape_rate function. 
This is useful for testing the function with different atmospheric compositions, and for making it easier to use the function in a more general context."""
import numpy as np

class Atmosphere:
    ### Dictionary of the most relevant dominant species for the wind, which comes with their own values for hte microphysics ###
    wind_microphysics = {
        #NB, these are the values for the singly ionised and fully dissociated equivalents!!!!!
        'H2': {'nu_0': 3.288467085473 * 10**15, 'mu_wind': 0.5, 'mu_plus_wind': 1},
        'H': {'nu_0': 3.288467085473 * 10**15, 'mu_wind': 0.5, 'mu_plus_wind': 1}, 
        'H2O': {'nu_0': 4.835981008048 * 10**15 , 'mu_wind': 3, 'mu_plus_wind': 6},
        'CO2' : {'nu_0' : 3.3368 * 10**15, 'mu_wind': 7.2, 'mu_plus_wind': 14.4}, 
        'N2' : {'nu_0' : 3.7721 * 10**15, 'mu_wind': 7, 'mu_plus_wind': 14}
    }

    def __init__(self, M_p, R_p, pressures, temperatures, heights, F_xuv=None, F_ins=None, dominant_species=None, T_wind=10**(4), vmrs=None, P_base=10**(-4), nu_0=None, mu_wind=None, mu_plus_wind=None):
        #input chosen by user
        self.P_base = P_base #[Pa] pressure at the base of the escaping atmosphere, REFERENCE Lopez et. al. 2017
        
        #from planet bulk proerties
        self.M_p = M_p
        self.R_p = R_p
        
        if F_xuv is not None:
            self.F_xuv = F_xuv #[W/m^2] XUV flux
        elif F_ins is not None:
            self.F_xuv = self.calc_xuv_from_insolation(F_ins) #[W/m^2] XUV flux calculated from total insolation
        else:
            raise ValueError("Either F_xuv or F_ins must be provided as an input to the Atmosphere class.")

        #standard wind properties
        self.T_wind = T_wind
       
        #imported atmospheric profiles of different properties as a function of radius
        self.T = temperatures
        self.radii = R_p + heights
        self.pressures = pressures
        self.vmrs = vmrs if vmrs is not None else None
        
        self.read_off_wind_base_parameters()    
        self.determine_wind_microphysics(nu_0, mu_wind, mu_plus_wind)
        

               

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
        self.P_base_at_R_base = self.pressures[index]
        
        self.R_base = self.radii[index] #[m] radius of the base of the escaping atmosphere, where radii: array of radii from the planetary radius to 10 times the planetary radius, index: index of minimum element from the array of the absolute difference between the pressure profile and the pressure at the base of the escaping atmosphere
        self.T_base = self.T[index] #[K] temperature at the base of the escaping atmosphere, where self.T: array of temperatures as a function of radius, index: index of the radius of the base of the escaping atmosphere
        
        if self.vmrs is not None:
            #find vmrs at base
            self.vmrs_base = {species: vmr_array[index] for species, vmr_array in self.vmrs.items()} #vmrs at the base of the escaping atmosphere, where self.vmrs: dictonary of arrays of volume mixing ratios as a function of radius, index: index of the radius of the base of the escaping atmosphere
            #find dominant species at base
            self.dominant_species = max(self.vmrs_base, key=self.vmrs_base.get)
            #print(f"Dominant species at the base of the escaping atmosphere: {self.dominant_species} with VMR of {self.vmrs_base[self.dominant_species]:.2e}")
        

    def determine_wind_microphysics(self, nu_0, mu_wind, mu_plus_wind):
        '''
        Determines the microphysics parameters for the escaping wind based on the dominant species at the base of the escaping atmosphere.

        Takes input parameters: 
        Atmosphere object, must have a vmrs attribute to determine the dominant species at the base of the escaping atmosphere from.

        All calculations done in SI units.
        '''
        #Checks for manual input
        if nu_0 is not None and mu_wind is not None and mu_plus_wind is not None:
            self.dominant_species =  "User_Defined"
            self.nu_0 = nu_0
            self.mu_wind = mu_wind
            self.mu_plus_wind = mu_plus_wind
            self.dominant_species_found_in_dict = False
            return
        elif self.vmrs is not None:
            #Otherwise use proteus data
            spec = self.dominant_species
            if spec in self.wind_microphysics:
                self.dominant_species_found_in_dict = True
                self.nu_0 = self.wind_microphysics[spec]['nu_0']
                self.mu_wind = self.wind_microphysics[spec]['mu_wind']
                self.mu_plus_wind = self.wind_microphysics[spec]['mu_plus_wind']
                return
            else:
                print(f'{self.dominant_species} found to be dominant, but not defined in microphysics dictionary. \n Add parameters (nu_0, mu_wind, mu_plus_wind) there or provide them manually.')
        else:
            raise ValueError("Microphysics parameters (nu_0, mu_wind, mu_plus_wind could not be determined. Please input vmrs or the values manually.)")

        

