from rr_escape import *

### Importing the atmospheres we want to examine ###

from proteus_fetch import atm_H2_case1
from atmospheres.simple_H2 import atm_H2


### Calculating the escape parameters ###

#I want to examine this atmosphere:
atm = atm_H2_case1
P_base = 10**(-4) #[Pa] pressure at the base of the escaping atmosphere, REFERENCE Lopez et. al. 2017

print(get_rr_escape_diagnostics(atm, P_base))