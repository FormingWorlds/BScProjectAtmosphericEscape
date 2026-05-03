import pandas as pd
import numpy as np
from scipy.constants import k, G

df = pd.read_csv("diffusion_coefficients.csv", comment="#")
#select a system, remember that this gives the diffusion coeff in cm^2/sec
row = df[df["system"] == "H_in_CO2"].iloc[0]
A = row["A"]
s = row["s"]

n_all = 100 #total number density, cm^-3
n_minor = 1 #number density of minor constituent, cm^-3
T = 300 #temperature profile, Kelvin
M = 1 #mass of the planet, kg
m_a = 1 #mean molecular mass, ?

T_0 = 1
r_0 = 1

mole_frac = n_minor/(n_all+n_minor) #mole fraction of the minor constituent
D12 = (A*(T**s))/n_all #diffusion coefficient, cm^2/sec
b = D12*n_all #binary diffusion param
H_a = k*T_0*(r_0**2)/G*M*m_a