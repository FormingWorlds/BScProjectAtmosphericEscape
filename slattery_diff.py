import numpy as np

#method from Reid page 271, "Slattery method"

molar_masses = { #[g/mol]
    "H2" : 2.016,  
    "H2O": 18.015,
    "H"  : 1.008,
    "N2": 28.014
}

critical_temps = {
    "H2O" : 647.1, #critial temp of water, [K]     NIST
    "H2" : 33.18, #critical temp of H2, [K]    NIST
    "H" : 33.18, #critical temp of H, [K]    NIST
    "N2" : 126.21 #[K], wiki
}

critical_pressures = {
    "H2O" : 217.75, #critical pressure, [atm]
    "H2" : 12.68986, #critical pressure,[atm] (1.2858 MPa, wiki)
    "H" : 12.68986, #critical pressure,[atm] (1.2858 MPa, wiki)
    "N2" : 33.4567 #[atm]
}

#first val for interdiffusion of nonpolar gases and for self-diffusion, second val for interdiffusion of water and nonpolar gas
k_c = np.array([2.74e-4, 3.63e-4])

#first val for interdiffusion of nonpolar gases and for selfdiffusion, second val for interdiffusion of water and nonpolar gas
n_c = np.array([1.823, 2.334])


def slattery_diff(T, P, minor, major):
    """calculates binary diffusion coefficient (D) for untabulated nonpolar-nonpolar (self also possible) or nonpolar-water species
    given temperature, pressure, all caps strings of major and minor constituents"""
    
    if ((minor == 'H2O') or (major == 'H2O')):
        k = k_c[1]
        n = n_c[1]
    else:
        k = k_c[0]
        n = n_c[0]
        
    Tc1 = critical_temps['minor']
    Tc2 = critical_temps['major']
    Pc1 = critical_pressures['minor']
    Pc2 = critical_pressures['major']
    M1 = molar_masses['minor']
    M2 = molar_masses['major']
    
    T_red = T / ((Tc1 * Tc2)**0.5)
    D = (k * (T_red**n) * ((Pc1 * Pc2)**0.33) * ((Tc1 * Tc2)**(5/12)) * ((M1 + M2)**0.5)) / (P * ((M1 * M2)**0.5))
    #D is in [cm^2/s]
    
    return (D)
    