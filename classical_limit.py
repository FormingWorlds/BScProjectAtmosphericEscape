import pandas as pd
import numpy as np

df = pd.read_csv("diffusion_coefficients.csv", comment="#")
#select a system, remember that this gives the diffusion coeff in cm^2/sec
row = df[df["system"] == "H_in_CO2"].iloc[0]
A = row["A"]
s = row["s"]