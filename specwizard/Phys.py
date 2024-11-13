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


  
    phys_dic = {
        "c": 29979245800.0,
        "Grav": 6.674079999999999e-08,
        "mH": 1.672621898e-24,
        "me": 9.10938356e-28,
        "amu": 1.66053904e-24,
        "Y": 0.248,
        "planck": 6.62607004e-27,
        "kB": 1.38064852e-16,
        "sigma": 6.63e-18,
        "sigmaT": 6.6524587158e-25,
        "alphaB": 3.334168827340852e-13,
        "eps0": 0.07957747154594767,
        "e": 4.80320467299766e-10,
        "GMsun": 1.3271244e+26,
        "Msun": 1.9884754153381438e+33,
        "pc": 3.0856775814913674e+18,
        "kpc": 3.0856775814913673e+21,
        "Mpc": 3.0856775814913676e+24,
        "micron": 0.0001,
        "yr": 31557600.0,
        "H100": 3.2407792894443648e-18,
        "eV": 1.6021766208000001e-12,
        "Ryd": 2.1798723253902544e-11,
        "ElementSolarMassFractions": {
            "info": "From Wiersma et al. 2009",
            "H": 0.7065,
            "He": 0.2806,
            "C": 2.07e-3,
            "N": 8.36e-4,
            "O": 5.49e-3,
            "Ne": 1.41e-3,
            "Mg": 5.91e-4,
            "Si": 6.83e-4,
            "S": 4.09e-4,
            "Ca": 6.44e-5,
            "Fe": 1.10e-3,
        },
    }
    return phys_dic
