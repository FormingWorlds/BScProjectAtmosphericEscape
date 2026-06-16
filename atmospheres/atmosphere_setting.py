"""This file makes a class that can hold atmosphere parameters such that they don't have to be input as separate arguments to the rr_escape_rate function. 
This is useful for testing the function with different atmospheric compositions, and for making it easier to use the function in a more general context."""
import numpy as np
import scipy.constants as sp
from radiation_recombination_coefficient_cauldron import species_rr_coefficients_case_B
from el_escape import evaporation_efficiencies

### Dictionary of the most relevant dominant species for the wind, which comes with their own values for hte microphysics ###
wind_microphysics = {
        #NB, these are the values for the singly ionised and fully dissociated equivalents!!!!!
        'H2': {'nu_0': 3.288467085473 * 10**15, 'mu_wind': 0.5, 'mu_plus_wind': 1, 'rr_coeff': species_rr_coefficients_case_B['H2'], 'epsilon_xuv': evaporation_efficiencies['H2']},
        'H': {'nu_0': 3.288467085473 * 10**15, 'mu_wind': 0.5, 'mu_plus_wind': 1, 'rr_coeff': species_rr_coefficients_case_B['H'], 'epsilon_xuv': evaporation_efficiencies['H']}, 
        'H2O': {'nu_0': 3.293303066481 * 10**15 , 'mu_wind': 3, 'mu_plus_wind': 6, 'rr_coeff': species_rr_coefficients_case_B['H2O'], 'epsilon_xuv': evaporation_efficiencies['H2O']},
        'CO2' : {'nu_0' : 3.293303066481 * 10**15, 'mu_wind': 7.333, 'mu_plus_wind': 14.667, 'rr_coeff': species_rr_coefficients_case_B['CO2'], 'epsilon_xuv': evaporation_efficiencies['CO2']}, 
        'N2' : {'nu_0' : 3.513340202347 * 10**15, 'mu_wind': 7, 'mu_plus_wind': 14, 'rr_coeff': species_rr_coefficients_case_B['N2'], 'epsilon_xuv': evaporation_efficiencies['N2']},
        'H2He' : {'nu_0' : 5.948256639899 * 10**15, 'mu_wind': 0.62, 'mu_plus_wind': 1.3 , 'rr_coeff': species_rr_coefficients_case_B['H2He'], 'epsilon_xuv': evaporation_efficiencies['H2He']},
    }
class Atmosphere:
    
    def __init__(self, M_p, pressures, heights, temperatures, F_xuv=None, dominant_species=None, T_wind=10**(4), vmrs=None, P_base=10**(-4), nu_0=None, mu_wind=None, mu_plus_wind=None, R_p=None, determine_radius=False, rr_coeff=None, epsilon_xuv=None):
        #input chosen by user
        self.P_base = P_base #[Pa] pressure at the base of the escaping atmosphere, REFERENCE Lopez et. al. 2017
        
        #from planet bulk proerties
        self.M_p = M_p
        
        if determine_radius is False and R_p is None:
            raise ValueError('Planetary radius not provided and option to calculate it based on mass and composition is False. Either set "determine_radius=True" or provide a planetary radius "R_p=value"')
        elif determine_radius is False and R_p is not None:
            self.R_p = R_p
        else: 
            self.R_p= self.determine_radius_from_MR_relation(self.M_p)


        if F_xuv is not None:
            self.F_xuv = F_xuv #[W/m^2] XUV flux
        else:
            raise ValueError("F_xuv  must be provided as an input to the Atmosphere class to calculate the escape rate.")

        #standard wind properties
        self.T_wind = T_wind
       
        #imported atmospheric profiles of different properties as a function of radius
        self.T = temperatures
        self.radii = R_p + heights
        self.pressures = pressures
        self.vmrs = vmrs if vmrs is not None else None
                
        self.read_off_wind_base_parameters()    
        self.determine_wind_microphysics(nu_0, mu_wind, mu_plus_wind, dominant_species, rr_coeff, epsilon_xuv)
        
        # First check if user provided a recombination coefficient manually, if not check if the dominant species is in the microphysics dictionary, else raise an error
        if rr_coeff is not None:
            self.rr_coeff = rr_coeff
        elif self.dominant_species_found_in_dict:
            self.rr_coeff = wind_microphysics[self.dominant_species]['rr_coeff']
        else:
            raise ValueError("No valid recombination coefficient found for the specified dominant species. Please provide a recombination coefficient manually or check that the dominant species is correctly specified and present in the microphysics dictionary.")


    
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
        

    def determine_wind_microphysics(self, nu_0, mu_wind, mu_plus_wind, dominant_species, rr_coeff, epsilon_xuv):
        '''
        Determines the microphysics parameters for the escaping wind based on the dominant species at the base of the escaping atmosphere.

        Takes input parameters: 
        Atmosphere object, must have a vmrs attribute to determine the dominant species at the base of the escaping atmosphere from. 

        All calculations done in SI units.
        '''
        #Checks for manual input first
        if nu_0 is not None and mu_wind is not None and mu_plus_wind is not None and rr_coeff is not None:
            if dominant_species is not None:
                self.dominant_species = dominant_species
            else: 
                self.dominant_species = 'User_defined'

            self.nu_0 = nu_0
            self.mu_wind = mu_wind
            self.mu_plus_wind = mu_plus_wind
            self.rr_coeff = rr_coeff
            self.epsilon_xuv = epsilon_xuv
            self.dominant_species_found_in_dict = False
            return
        
        # If there is no manual input, we check if this is proteus data with vmrs which can give us the dominant species and microphysics parameters
        elif self.vmrs is not None:
            #Otherwise use proteus data
            spec = self.dominant_species
            if spec in wind_microphysics:
                self.dominant_species_found_in_dict = True
                self.nu_0 = wind_microphysics[spec]['nu_0']
                self.mu_wind = wind_microphysics[spec]['mu_wind']
                self.mu_plus_wind = wind_microphysics[spec]['mu_plus_wind']
                self.rr_coeff = wind_microphysics[spec]['rr_coeff']
                self.epsilon_xuv = wind_microphysics[spec]['epsilon_xuv']
                return
            else:
                raise ValueError(f'{self.dominant_species} found to be dominant, but not defined in microphysics dictionary. \n Add parameters (nu_0, mu_wind, mu_plus_wind, rr_coeff) there or provide them manually.')

        # And if there is no manual microphysics and this isn't a proteus atmosphere, we check is there is input on the dominant species to try to assign microphysics parameters that way
        elif dominant_species is not None:
            if dominant_species in wind_microphysics:
                self.dominant_species_found_in_dict = True
                self.dominant_species = dominant_species
                self.nu_0 = wind_microphysics[dominant_species]['nu_0']
                self.mu_wind = wind_microphysics[dominant_species]['mu_wind']
                self.mu_plus_wind = wind_microphysics[dominant_species]['mu_plus_wind']
                self.rr_coeff = wind_microphysics[dominant_species]['rr_coeff']
                self.epsilon_xuv = wind_microphysics[dominant_species]['epsilon_xuv']
            else:
                raise ValueError(f'{dominant_species} specified as dominant, but not defined in microphysics dictionary. \n Add parameters (nu_0, mu_wind, mu_plus_wind, rr_coeff) there or provide them manually.')
        else:
            raise ValueError("Microphysics parameters (nu_0, mu_wind, mu_plus_wind, rr_coeff, epsilon_xuv) could not be determined. Please input vmrs or the values manually.)")

    @staticmethod
    def determine_radius_from_MR_relation(M_p):
        M_earth = 5.9722 * 10**24  # [kg]
        R_earth = 6.371 * 10**6    # [m] (Mean Earth radius in SI meters)
    
        #Convert input mass from kg to Earth masses
        M_p_earth = M_p / M_earth
    
        ### Mass-radius relation adapted from Parc et. al. 2024 ###
        if M_p_earth < 10:
            R_p_earth = 1.02 * M_p_earth**0.28
        
        elif 10 <= M_p_earth < 138:
            R_p_earth = 0.61 * M_p_earth**0.67
        
        else:  # aka M_p_earth >= 138
            R_p_earth = 11.9 * M_p_earth**0.01
        
        #Convert the final radius from Earth radii back to SI meters and return
        return R_p_earth * R_earth

    def determine_temperature_from_bolometric_flux(F_bol):
        '''
        Determines the effective temperature of the planet based on the bolometric flux and planetary radius from the Stefan-Boltzmann law and assuming the planet is a black body.

        Takes input parameters: F_bol [W/m^2] - bolometric flux, R_p [m] - planetary radius

        Output: T_eff [K] - effective temperature of the planet
        '''
        sigma_sb = sp.constants.sigma #[W m^-2 K^-4] Stefan-Boltzmann constant
        T_eq = (F_bol / (4 * sigma_sb))**(1/4) #[K] effective temperature of the planet, where F_bol: bolometric flux, R_p: planetary radius, sigma_sb: Stefan-Boltzmann constant
        return T_eq

