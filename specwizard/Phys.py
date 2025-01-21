import unyt
def ReadPhys():
    """
    Reads and returns a dictionary of physical constants and parameters commonly used in astrophysical calculations,
    including fundamental constants, solar mass fractions, and units in CGS (centimeter-gram-second) system.

    Returns:
        dict: A dictionary containing key-value pairs for various physical constants:
            - Fundamental constants:
                - "c" (float): Speed of light in cm/s.
                - "Grav" (float): Gravitational constant in cm^3 g^-1 s^-2.
                - "mH" (float): Mass of a hydrogen atom in g.
                - "me" (float): Mass of an electron in g.
                - "amu" (float): Atomic mass unit in g.
                - "planck" (float): Planck constant in erg s.
                - "kB" (float): Boltzmann constant in erg K^-1.
                - "sigma" (float): Stefan-Boltzmann constant in erg cm^-2 s^-1 K^-4.
                - "sigmaT" (float): Thomson cross-section in cm^2.
                - "alphaB" (float): Recombination coefficient in cm^3 s^-1.
                - "eps0" (float): Permittivity of free space in cm^-3 g s^2 C^-2.
                - "e" (float): Elementary charge in statcoulombs (esu).
            - Solar and galactic constants:
                - "GMsun" (float): Gravitational constant times solar mass in cm^3 s^-2.
                - "Msun" (float): Solar mass in g.
                - "pc" (float): Parsec in cm.
                - "kpc" (float): Kiloparsec in cm.
                - "Mpc" (float): Megaparsec in cm.
                - "micron" (float): Micron in cm.
                - "yr" (float): Year in seconds.
                - "H100" (float): Hubble constant in s^-1 (for H0 = 100 km/s/Mpc).
                - "eV" (float): Electronvolt in erg.
                - "Ryd" (float): Rydberg constant in erg.
            - Elemental solar mass fractions (from Wiersma et al. 2009):
                - "ElementSolarMassFractions" (dict): Dictionary containing:
                    - "H" (float): Hydrogen mass fraction.
                    - "He" (float): Helium mass fraction.
                    - "C" (float): Carbon mass fraction.
                    - "N" (float): Nitrogen mass fraction.
                    - "O" (float): Oxygen mass fraction.
                    - "Ne" (float): Neon mass fraction.
                    - "Mg" (float): Magnesium mass fraction.
                    - "Si" (float): Silicon mass fraction.
                    - "S" (float): Sulfur mass fraction.
                    - "Ca" (float): Calcium mass fraction.
                    - "Fe" (float): Iron mass fraction.

    Example:
        >>> constants = ReadPhys()
        >>> print(constants["c"])
        29979245800.0
    """


    length_cgs = unyt.cm
    time_cgs   = unyt.s
    velocity_cgs = length_cgs/time_cgs
    mass_cgs     = unyt.g
    dimmensionless = unyt.dimensionless
    energy_cgs   = mass_cgs * (velocity_cgs)**2
    temperature_cgs = unyt.K

    phys_dic = {
        "c": 29979245800.0*velocity_cgs,
        "Grav": 6.674079999999999e-08 *((length_cgs)**3 * (mass_cgs)**-1* (time_cgs)**-2),
        "mH": 1.672621898e-24*mass_cgs,
        "me": 9.10938356e-28*mass_cgs,
        "amu": 1.66053904e-24*mass_cgs,
        "Y": 0.248*dimmensionless,
        "planck": 6.62607004e-27 * energy_cgs*time_cgs,
        "kB": 1.38064852e-16 * energy_cgs/temperature_cgs,
        "sigma": 6.63e-18 * energy_cgs * length_cgs**-2 * time_cgs**-1 * temperature_cgs**-4,
        "sigmaT": 6.6524587158e-25 * length_cgs**2,
        "alphaB": 3.334168827340852e-13 * length_cgs**3* time_cgs**-1,
        "eps0": 0.07957747154594767,
        "e": 4.80320467299766e-10*unyt.statC,
        "GMsun": 1.3271244e+26,
        "Msun": 1.9884754153381438e+33*mass_cgs,
        "pc": 3.0856775814913674e+18*length_cgs,
        "kpc": 3.0856775814913673e+21*length_cgs,
        "Mpc": 3.0856775814913676e+24*length_cgs,
        "micron": 0.0001*length_cgs,
        "yr": 31557600.0*time_cgs,
        "H100": 3.2407792894443648e-18*time_cgs**-1,
        "eV": 1.6021766208000001e-12*energy_cgs,
        "Ryd": 2.1798723253902544e-11*energy_cgs,
        "ElementSolarMassFractions": {
            "info": "From Wiersma et al. 2009",
            "H": 0.7065*dimmensionless,
            "He": 0.2806*dimmensionless,
            "C": 2.07e-3*dimmensionless,
            "N": 8.36e-4*dimmensionless,
            "O": 5.49e-3*dimmensionless,
            "Ne": 1.41e-3*dimmensionless,
            "Mg": 5.91e-4*dimmensionless,
            "Si": 6.83e-4*dimmensionless,
            "S": 4.09e-4*dimmensionless,
            "Ca": 6.44e-5*dimmensionless,
            "Fe": 1.10e-3*dimmensionless,
        },
    }
    return phys_dic
