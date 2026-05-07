### Importing the escape functions ###
from rr_escape import *

### Importing the single atmospheres we want to examine ###
from proteus_fetch import atm_H2_case1
#from proteus_fetch import atm_H2_case1000
#from atmospheres.simple_H2 import atm_H2
#from atmospheres.simple_H2He_mix import atm_H2He_mix
#from atmospheres.simple_H2O import atm_H2O


### Doing escape calculations ###
atm = atm_H2_case1
examine_atmosphere_for_rr_escape(atm)
