### Importing the escape functions ###
from rr_escape import *

### Importing the single atmospheres we want to examine ###
from proteus_fetch import proteus_atmosphere_cases

from atmospheres.simple_H2 import atm_H2
from atmospheres.simple_H2He_mix import atm_H2He_mix
from atmospheres.simple_H2O import atm_H2O


### Doing escape calculations ###
#print(dict.keys(proteus_atmosphere_cases))
#atm = proteus_atmosphere_cases['H2_1_M_earth_1_F_earth']

examine_atmosphere_for_rr_escape(atm_H2)
examine_atmosphere_for_rr_escape(atm_H2He_mix)
examine_atmosphere_for_rr_escape(atm_H2O)

#atm_H2_1M_1000F = proteus_atmosphere_cases['H2_1_M_earth_1000_F_earth']
#atm_H2_10M_1F = proteus_atmosphere_cases['H2_10_M_earth_1_F_earth']
#atm_H2O_1M_1000F = proteus_atmosphere_cases['H2O_1_M_earth_1000_F_earth']

#examine_atmosphere_for_rr_escape(atm_H2_1M_1000F)
#examine_atmosphere_for_rr_escape(atm_H2O_1M_1000F)




