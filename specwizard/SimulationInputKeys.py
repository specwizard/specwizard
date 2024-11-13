def get_simkeys(sim_name='eagle'):
    """
    Retrieves a dictionary of key mappings for the specified simulation type, which includes element names, group names,
    and various property labels used in different astrophysical simulations such as 'eagle', 'swift', 'hydrangea', 'colibre',
    and 'illustris'.

    Args:
        sim_name (str): The name of the simulation type. Options are:
                        - "eagle"
                        - "swift"
                        - "hydrangea"
                        - "colibre"
                        - "illustris"
                        Default is 'eagle'.

    Returns:
        dict: A dictionary of keys specific to the simulation type with the following structure:
            - `elementnames` (list): List of element names.
            - `groupname` (str): Name of the particle or line-of-sight group.
            - Various properties (str): Labels for attributes such as densities, masses, positions, velocities, temperatures,
              metallicities, star formation rate, ion fractions, etc.

    Simulation Types:
        - eagle: Keys for 'eagle' simulation data with `elementnames`, 'Header', and 'los' and 'snapshot' keys.
        - swift: Keys for 'swift' simulation data with `elementnames`, 'Header', 'los', and 'snapshot' keys.
        - hydrangea: Keys for 'hydrangea' simulation data with `elementnames`, 'Header', and 'snapshot' keys.
        - colibre: Keys for 'colibre' simulation data with `elementnames`, 'Header', 'los', and 'snapshot' keys.
        - illustris: Keys for 'illustris' simulation data with `elementnames` and 'snapshot' keys only.

    Example:
        >>> keys = get_simkeys('eagle')
        >>> print(keys['los']['elementnames'])
        ['Hydrogen', 'Helium', 'Nitrogen', 'Carbon', 'Oxygen', 'Magnesium', 'Silicon', 'Iron']

    Raises:
        KeyError: If an unsupported simulation type is passed as `sim_name`.
    """
    
    elementnames_eagle = ['Hydrogen', 'Helium', 'Nitrogen', 'Carbon', 'Oxygen', 'Magnesium', 'Silicon', 'Iron']
    elementnames_hydra = ['Hydrogen', 'Helium', 'Nitrogen', 'Carbon', 'Oxygen', 'Magnesium', 'Silicon', 'Iron']
    elementnames_swift = ['Hydrogen', 'Helium', 'Carbon', 'Nitrogen', 'Oxygen', 'Neon', 'Magnesium', 'Silicon', 'Iron', 'Europium']
    elementnames_TNG   = ['Hydrogen', 'Helium', 'Carbon', 'Nitrogen', 'Oxygen', 'Neon', 'Magnesium', 'Silicon', 'Iron']
    
    sym_keys = {
        'eagle': {
            'Header': {'Redshift': ['Header', 'Redshift']},
            'los': {
                'elementnames': elementnames_eagle,
                'groupname': 'LOS{}',
                'ElementAbundance': 'ElementAbundance',
                'Densities': 'Density',
                'SmoothingLengths': 'SmoothingLength',
                'Masses': 'Mass',
                'Positions': 'Positions',
                'Velocities': 'Velocity',
                'Temperatures': 'Temperature',
                'Metallicities': 'Metallicity',
                'StarFormationRate': 'StarFormationRate',
                'x-axis': 'x-axis',
                'y-axis': 'y-axis',
                'z-axis': 'z-axis',
                'x-position': 'x-position',
                'y-position': 'y-position',
                'z-position': 'z-position',
                'Number_of_part_this_los': 'Number_of_part_this_los'
            },
            'snapshot': {
                'elementnames': elementnames_eagle,
                'groupname': 'PartType0',
                'ElementAbundance': 'ElementAbundance',
                'Densities': 'Density',
                'SmoothingLengths': 'SmoothingLength',
                'Masses': 'Mass',
                'Positions': 'Coordinates',
                'Velocities': 'Velocity',
                'Temperatures': 'Temperature',
                'Metallicities': 'Metallicity',
                'StarFormationRate': 'StarFormationRate'
            }
        },
        'swift': {
            'Header': {'Redshift': ['Cosmology', 'Redshift']},
            'los': {
                'elementnames': elementnames_swift,
                'groupname': 'LOS_{:04d}',
                'ElementAbundance': 'ElementMassFractions',
                'Densities': 'Densities',
                'SmoothingLengths': 'SmoothingLengths',
                'Masses': 'Masses',
                'Positions': 'Coordinates',
                'Velocities': 'Velocities',
                'Temperatures': 'Temperatures',
                'Metallicities': 'MetalMassFractions',
                'StarFormationRate': 'StarFormationRates',
                'x-axis': 'Xaxis',
                'y-axis': 'Yaxis',
                'z-axis': 'Zaxis',
                'x-position': 'Xpos',
                'y-position': 'Ypos',
                'z-position': 'Ypos',
                'Number_of_part_this_los': 'NumParts'
            },
            'snapshot': {
                'elementnames': elementnames_swift,
                'groupname': 'PartType0',
                'ElementAbundance': 'ElementMassFractions',
                'Densities': 'Densities',
                'SmoothingLengths': 'SmoothingLengths',
                'Masses': 'Masses',
                'Positions': 'Coordinates',
                'Velocities': 'Velocities',
                'Temperatures': 'Temperatures',
                'Metallicities': 'MetalMassFractions',
                'StarFormationRate': 'StarFormationRates'
            }
        },
        'hydrangea': {
            'Header': {'Redshift': ['Header', 'Redshift']},
            'snapshot': {
                'elementnames': elementnames_hydra,
                'groupname': 'PartType0',
                'ElementAbundance': 'ElementAbundance',
                'Densities': 'Density',
                'SmoothingLengths': 'SmoothingLength',
                'Masses': 'Mass',
                'Positions': 'Coordinates',
                'Velocities': 'Velocity',
                'Temperatures': 'Temperature',
                'Metallicities': 'Metallicity',
                'StarFormationRate': 'StarFormationRate'
            }
        },
        'colibre': {
            'Header': {'Redshift': ['Cosmology', 'Redshift']},
            'los': {
                'elementnames': elementnames_swift,
                'groupname': 'LOS_{:04d}',
                'ElementAbundance': 'ElementMassFractionsInGas',
                'Densities': 'Densities',
                'SmoothingLengths': 'SmoothingLengths',
                'Masses': 'Masses',
                'Positions': 'Coordinates',
                'Velocities': 'Velocities',
                'Temperatures': 'Temperatures',
                'Metallicities': 'MetalMassFractions',
                'StarFormationRate': 'StarFormationRates',
                'IonFractions': 'SpeciesFractions',
                'x-axis': 'Xaxis',
                'y-axis': 'Yaxis',
                'z-axis': 'Zaxis',
                'x-position': 'Xpos',
                'y-position': 'Ypos',
                'z-position': 'Ypos',
                'Number_of_part_this_los': 'NumParts'
            },
            'snapshot': {
                'elementnames': elementnames_swift,
                'groupname': 'PartType0',
                'ElementAbundance': 'ElementMassFractions',
                'Densities': 'Densities',
                'SmoothingLengths': 'SmoothingLengths',
                'Masses': 'Masses',
                'Positions': 'Coordinates',
                'Velocities': 'Velocities',
                'Temperatures': 'Temperatures',
                'Metallicities': 'MetalMassFractions',
                'StarFormationRate': 'StarFormationRates',
                'IonFractions': 'SpeciesFractions'
            }
        },
        'illustris': {
            'snapshot': {
                'elementnames': elementnames_TNG,
                'groupname': 'PartType0',
                'ElementAbundance': 'GFM_Metals',
                'Densities': 'Density',
                'Masses': 'Masses',
                'Positions': 'Coordinates',
                'Velocities': 'Velocities',
                'InternalEnergy': 'InternalEnergy',
                'ElectronAbundance': 'ElectronAbundance',
                'Metallicities': 'GFM_Metallicity',
                'IonFractions': 'SpeciesFractions',
                'NeutralHydrogenAbundance': 'NeutralHydrogenAbundance'
            }
        }
    }

    return sym_keys[sim_name]

    