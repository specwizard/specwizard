
def get_simkeys(sim_name='eagle'):
                
    elementnames_eagle = ['Hydrogen', 'Helium', 'Nitrogen', 'Carbon', 'Oxygen', 'Magnesium', 'Silicon','Iron']
    elementnames_hydra = ['Hydrogen', 'Helium', 'Nitrogen', 'Carbon', 'Oxygen', 'Magnesium', 'Silicon','Iron']
    elementnames_swift = ['Hydrogen', 'Helium','Carbon','Nitrogen','Oxygen','Neon','Magnesium','Silicon','Iron','Europium']
    elementnames_TNG   = ['Hydrogen', 'Helium','Carbon','Nitrogen','Oxygen','Neon','Magnesium','Silicon','Iron']
    sym_keys = {}
    sym_keys['eagle'] = {}
    sym_keys['eagle']['Header'] = {'Redshift': ['Header','Redshift']}
    sym_keys['eagle']['los'] ={'elementnames'           : elementnames_eagle,
                               'groupname'              : 'LOS{}',
                               'ElementAbundance'       : 'ElementAbundance',
                               'Densities'              : 'Density',
                               'SmoothingLengths'       : 'SmoothingLength',
                               'Masses'                 : 'Mass',
                               'Positions'              : 'Positions',
                               'Velocities'             : 'Velocity',
                               'Temperatures'           : 'Temperature',
                               'Metallicities'          : 'Metallicity',
                               'StarFormationRate'      : 'StarFormationRate',
                               'x-axis'                 : 'x-axis',
                               'y-axis'                 : 'y-axis',
                               'z-axis'                 : 'z-axis',
                               'x-position'             : 'x-position',
                               'y-position'             : 'y-position',
                               'z-position'             : 'z-position',
                               'Number_of_part_this_los': 'Number_of_part_this_los'}


    sym_keys['eagle']['snapshot'] =  { 'elementnames'      : elementnames_eagle,
                                       'groupname'         : 'PartType0',
                                       'ElementAbundance'  : 'ElementAbundance',
                                       'Densities'         : 'Density',
                                       'SmoothingLengths'  : 'SmoothingLength',
                                       'Masses'            : 'Mass',
                                       'Positions'         : 'Coordinates',
                                       'Velocities'        : 'Velocity',
                                       'Temperatures'      : 'Temperature',
                                       'Metallicities'     : 'Metallicity',
                                       'StarFormationRate' : 'StarFormationRate'}
    sym_keys['swift'] = {}
    sym_keys['swift']['Header'] = {'Redshift': ['Cosmology','Redshift']}

    sym_keys['swift']['los'] = {'elementnames'          : elementnames_swift,
                               'groupname'              : 'LOS_{:04d}',
                               'ElementAbundance'       : 'ElementMassFractions',
                               'Densities'              : 'Densities',
                               'SmoothingLengths'       : 'SmoothingLengths',
                               'Masses'                 : 'Masses',
                               'Positions'              : 'Coordinates',
                               'Velocities'             : 'Velocities',
                               'Temperatures'           : 'Temperatures',
                               'Metallicities'          : 'MetalMassFractions',
                               'StarFormationRate'      : 'StarFormationRates',
                               'x-axis'                 : 'Xaxis',
                               'y-axis'                 : 'Yaxis',
                               'z-axis'                 : 'Zaxis',
                               'x-position'             : 'Xpos',
                               'y-position'             : 'Ypos',
                               'z-position'             : 'Ypos',
                               'Number_of_part_this_los': 'NumParts'}


    sym_keys['swift']['snapshot'] = {'elementnames'     : elementnames_swift,
                               'groupname'              : 'PartType0',
                               'ElementAbundance'       : 'ElementMassFractions',
                               'Densities'              : 'Densities',
                               'SmoothingLengths'       : 'SmoothingLengths',
                               'Masses'                 : 'Masses',
                               'Positions'              : 'Coordinates',
                               'Velocities'             : 'Velocities',
                               'Temperatures'           : 'Temperatures',
                               'Metallicities'          : 'MetalMassFractions',
                               'StarFormationRate'      : 'StarFormationRates'}
    
    sym_keys['hydrangea'] = {}
    sym_keys['hydrangea']['Header'] = {'Redshift': ['Header','Redshift']}

    sym_keys['hydrangea']['snapshot'] =  {'elementnames'   : elementnames_hydra,
                                       'groupname'         : 'PartType0',
                                       'ElementAbundance'  : 'ElementAbundance',
                                       'Densities'         : 'Density',
                                       'SmoothingLengths'  : 'SmoothingLength',
                                       'Masses'            : 'Mass',
                                       'Positions'         : 'Coordinates',
                                       'Velocities'        : 'Velocity',
                                       'Temperatures'      : 'Temperature',
                                       'Metallicities'     : 'Metallicity',
                                       'StarFormationRate' : 'StarFormationRate'}
    sym_keys['colibre'] = {}
    sym_keys['colibre']['Header'] = {'Redshift': ['Cosmology','Redshift']}

    sym_keys['colibre']['los'] = {'elementnames'          : elementnames_swift,
                               'groupname'              : 'LOS_{:04d}',
                               'ElementAbundance'       : 'ElementMassFractionsInGas',
                               'Densities'              : 'Densities',
                               'SmoothingLengths'       : 'SmoothingLengths',
                               'Masses'                 : 'Masses',
                               'Positions'              : 'Coordinates',
                               'Velocities'             : 'Velocities',
                               'Temperatures'           : 'Temperatures',
                               'Metallicities'          : 'MetalMassFractions',
                               'StarFormationRate'      : 'StarFormationRates',
                               'x-axis'                 : 'Xaxis',
                               'y-axis'                 : 'Yaxis',
                               'z-axis'                 : 'Zaxis',
                               'x-position'             : 'Xpos',
                               'y-position'             : 'Ypos',
                               'z-position'             : 'Ypos',
                               'Number_of_part_this_los': 'NumParts',
                               'IonFractions'           : 'SpeciesFractions'}


    sym_keys['colibre']['snapshot'] = {'elementnames'     : elementnames_swift,
                               'groupname'              : 'PartType0',
                               'ElementAbundance'       : 'ElementMassFractions',
                               'Densities'              : 'Densities',
                               'SmoothingLengths'       : 'SmoothingLengths',
                               'Masses'                 : 'Masses',
                               'Positions'              : 'Coordinates',
                               'Velocities'             : 'Velocities',
                               'Temperatures'           : 'Temperatures',
                               'Metallicities'          : 'MetalMassFractions',
                               'StarFormationRate'      : 'StarFormationRates',
                               'IonFractions'           : 'SpeciesFractions'}
    

    sym_keys['illustris'] = {}
    sym_keys['illustris']['snapshot'] = {'elementnames' : elementnames_TNG,
                               'groupname'              : 'PartType0',
                               'ElementAbundance'       : 'GFM_Metals',
                               'Densities'              : 'Density',
                               'Masses'                 : 'Masses',
                               'Positions'              : 'Coordinates',
                               'Velocities'             : 'Velocities',
                               'InternalEnergy'         : 'InternalEnergy',
                               'ElectronAbundance'      : 'ElectronAbundance',
                               'Metallicities'          : 'GFM_Metallicity',
                               'IonFractions'           : 'SpeciesFractions',
                               'NeutralHydrogenAbundance': 'NeutralHydrogenAbundance'}


    return sym_keys[sim_name]

    