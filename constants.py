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

PARTICLE_RADII = {
    "H": 120e-12, #approximate
    "H2": 144.53e-12,
    "He": 130e-12,
    "C": 170e-12, #van der waals radius instead
    "N": 155e-12, #van der waals radius instead
    "O": 152e-12, #van der waals radius instead
    "N2": 182e-12,
    "O2": 173e-12,
    "CO": 188e-12,
    "CO2": 165e-12,
    "H2O": 132.5e-12 ,
    "OH": 150e-12, #approximate
    "NO": 158.5e-12,
    "NO2": 170e-12, #approximate
    "N2O": 165e-12,
    "CH4": 190e-12,
    "NH3": 130e-12,
    "H2S": 180e-12,
    "S": 180e-12, #van der waals radius instead
    "S2": 190e-12, #approximate
    "SO2": 180e-12,   

}