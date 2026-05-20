### Importing the escape functions ###
from rr_escape import examine_atmosphere_for_rr_escape, get_rr_escape_diagnostics

### Importing the single atmospheres we want to examine ###
from proteus_fetch import proteus_atmosphere_cases

from atmospheres.simple_H2 import atm_H2
from atmospheres.simple_H2He_mix import atm_H2He_mix
from atmospheres.simple_H2O import atm_H2O
import os
width = os.get_terminal_size().columns 


### Doing escape calculations ###
#print(dict.keys(proteus_atmosphere_cases))
#atm = proteus_atmosphere_cases['H2_1_M_earth_1_F_earth']
print('-' * width)
print('SIMPLE ATMOSPHERES')
print('-' * width)
print('H2')
examine_atmosphere_for_rr_escape(atm_H2)
print('-' * width)
print('H2He')
examine_atmosphere_for_rr_escape(atm_H2He_mix)
print('-' * width)
print('H2O')
examine_atmosphere_for_rr_escape(atm_H2O)
print('-' * width)


atm_H2_1M_1000F = proteus_atmosphere_cases['H2_1_M_earth_1000_F_earth']
#atm_H2_10M_1F = proteus_atmosphere_cases['H2_10_M_earth_1_F_earth']
atm_H2O_1M_1000F = proteus_atmosphere_cases['H2O_1_M_earth_1000_F_earth']
atm_CO2_1M_1000F = proteus_atmosphere_cases['CO2_1_M_earth_1000_F_earth']
atm_N2_1M_1000F = proteus_atmosphere_cases['N2_1_M_earth_1000_F_earth']

print('PROTEUS DATA')
print('-' * width)
print('H2')
examine_atmosphere_for_rr_escape(atm_H2_1M_1000F)
print('-' * width)
print('H2O')
examine_atmosphere_for_rr_escape(atm_H2O_1M_1000F)
print('-' * width)
print('CO2')
examine_atmosphere_for_rr_escape(atm_CO2_1M_1000F)
print('-' * width)
print('N2')
examine_atmosphere_for_rr_escape(atm_N2_1M_1000F)




