from scipy.constants import atomic_mass

R_EARTH = 6.371e6
M_EARTH = 5.972e24

SPECIES_MASSES = {
    "H": 1 * atomic_mass,
    "H2": 2 * atomic_mass,
    "He": 4 * atomic_mass,
    "C": 12 * atomic_mass,
    "N": 14 * atomic_mass,
    "O": 16 * atomic_mass,
    "N2": 28 * atomic_mass,
    "O2": 32 * atomic_mass,
    "CO": 28 * atomic_mass,
    "CO2": 44 * atomic_mass,
    "H2O": 18 * atomic_mass,
    "OH": 17 * atomic_mass,
    "NO": 30 * atomic_mass,
    "NO2": 46 * atomic_mass,
    "N2O": 44 * atomic_mass,
    "CH4": 16 * atomic_mass,
    "NH3": 17 * atomic_mass,
    "H2S": 34 * atomic_mass,
    "S": 32 * atomic_mass,
    "S2": 64 * atomic_mass,
    "SO2": 64 * atomic_mass,   
}