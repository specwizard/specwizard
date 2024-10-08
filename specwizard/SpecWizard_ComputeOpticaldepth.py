import numpy as np
from .SpecWizard_Elements import Elements
from .SpecWizard_Lines import Lines
from .Phys import ReadPhys
constants = ReadPhys()
#
import scipy.interpolate as interpolate

class ComputeOpticaldepth:
    ''' 
    Methods to compute optical depth as a function of velocity for a single sight line.
    
    Args:
        sightlineprojection (dict): A dictionary containing parameters related to the sight line projection.
    '''
    def __init__(self, sightlineprojection):
        '''
        Initializes the ComputeOpticaldepth class with the provided sight line projection.

        Args:
            sightlineprojection (dict): The specwizard dictionary containing parameters for the sight line projection,
                including extra parameters, element parameters, ion parameters, and optical depth parameters.

        Attributes:
            specparams (dict): Contains all the parameters from the sight line projection.
            periodic (bool): A flag indicating if the computation should be done with periodic boundary conditions.
            elements (list): List of element names from the element parameters.
            transitions (dict): Dictionary containing transition parameters for each ion.
            constants (dict): Physical constants read using the Phys.ReadPhys() method.
            ThermEff (bool): Indicates whether thermal effects are disabled (True if they are off).
            PecVelEff (bool): Indicates whether peculiar velocity effects are disabled (True if they are off).
            VoigtOff (bool): Indicates whether Voigt profile effects are disabled (True if they are off).
        '''
        self.specparams = sightlineprojection
        self.periodic   = self.specparams['extra_parameters']['periodic']
        
        # For each of the ions, determine rest-wavelength and f-value of the first transition
        self.elements = self.specparams["elementparams"]["ElementNames"]

        # Transition parameters for ions
        self.transitions = self.specparams["ionparams"]["transitionparams"]
        
        # Read physical constants
        self.constants  = Phys.ReadPhys()
        
        # Flags for various effects physical effects of the optical depth. 
        self.ThermEff   = self.specparams['ODParams']['ThermalEffectsOff']
        self.PecVelEff  = self.specparams['ODParams']['PecVelEffectsOff']
        self.VoigtOff   = self.specparams['ODParams']['VoigtOff']

    def MakeAllOpticaldepth(self, projected_los):
        ''' 
        Applies MakeOpticaldepth to compute optical depths for all desired ionic transitions.

        Args:
            projected_los (dict): The projection dictionary containing output from a sightline projection.
                - projected_los['SightInfo'] is obtained from snapshot.ReadParticles(sightline=sightline)["SightInfo"].
                - projected_los['Projection'] is an instance of SightLineProjection.

        Returns:
            dict: A dictionary containing the computed optical depth spectra for the various ionic transitions.

        Notes:
            This method computes the optical depths for the desired ionic transitions, taking into account different effects 
            such as thermal effects, peculiar velocities, and Voigt profile effects. 
        '''
        
        projection = {}
        projection["SightInfo"] = self.specparams["sightline"]
        projection["Header"] = self.specparams["Header"]
        projection["Projection"] = projected_los
        
        # Header information from the snapshot
        self.header = projection["Header"]
        
        # Pixel properties
        pixel_kms = self.ToCGS(self.header, projection["Projection"]["pixel_kms"]) / 1e5  # One pixel in km/s
        pixel = self.ToCGS(self.header, projection["Projection"]["pixel"])  # Size of pixel in cm
        sight_kms = self.ToCGS(self.header, projection["SightInfo"]["sightkms"]) / 1e5
        projection["SightInfo"]["Sight_kms"] = sight_kms
        
        self.sightinfo = projection["SightInfo"]
        
        # Velocity in km/s for each pixel
        vel_kms = np.arange(projection["Projection"]["npix"]) * pixel_kms
        vunit = self.SetUnit(vardescription='Hubble velocity', Lunit=1e5, aFact=0, hFact=0)
        npix = len(vel_kms)
        

        sightparams = {
            'sight_kms': sight_kms,
            'vel_kms': vel_kms,
            'pixel_kms': pixel_kms,
            'pixel': pixel
        }
        
        Ions = self.specparams["ionparams"]["Ions"]
        projectionIW = projection["Projection"]["Ion-weighted"]
        
        #If extend is true it will introduce the projected spectra into a larger zero's array so we can observe large effects of for exmaple a galaxy
        extend = self.specparams['sightline']['ProjectionExtend']
        if extend["extend"]:
            extended_factor = extend["extendfactor"]
            extended_npix = npix * extended_factor
            extended_vel_kms = np.arange(extended_npix) * pixel_kms
            start_indx = int(0.5 * npix * (extended_factor - 1))
            
            for ion in projectionIW.keys():
                for key in ['Densities', 'Velocities', 'Temperatures']:
                    temp_array = np.zeros_like(extended_vel_kms)
                    temp_array[start_indx:start_indx+npix] = projectionIW[ion][key]['Value'].copy()
                    projectionIW[ion][key]['Value'] = temp_array

            if 'SimIon-weighted' in projection["Projection"].keys():
                projectionSimIon = projection["Projection"]['SimIon-weighted']
                for ion in projectionSimIon.keys():
                    for key in ['Densities', 'Velocities', 'Temperatures']:
                        temp_array = np.zeros_like(extended_vel_kms)
                        temp_array[start_indx:start_indx+npix] = projectionSimIon[ion][key]['Value'].copy()
                        projectionSimIon[ion][key]['Value'] = temp_array

            sightparams['vel_kms'] = extended_vel_kms
            sightparams['sight_kms'] = extended_vel_kms.max()
        
        # Generate spectra
        spectra = self.WrapSpectra(Ions, projectionIW, sightparams, vel_mod=False, therm_mod=False)
        
        # Apply various effects
        if self.ThermEff:
            spectra['ThermalEffectsOff'] = self.WrapSpectra(Ions, projectionIW, sightparams, vel_mod=False, therm_mod=True)

        if self.PecVelEff:
            spectra['PeculiarVelocitiesOff'] = self.WrapSpectra(Ions, projectionIW, sightparams, vel_mod=True, therm_mod=False)

        if self.ThermEff and self.PecVelEff:
            spectra['ThermalPecVelOff'] = self.WrapSpectra(Ions, projectionIW, sightparams, vel_mod=True, therm_mod=True)

        DoSimIons = False
        #We will check if there is ion fraction calculations directly from the simulation and not from ion tables. 
        try:
            projectionSIW = projection["Projection"]['SimIon-weighted']
            if extend["extend"]:
                projectionSIW = projectionSimIon
            SimIons = list(projectionSIW.keys())
            DoSimIons = True
        except:
            pass

        if DoSimIons:
            SimIons = np.array(SimIons)
            all_ions = np.array(Ions)[:, 1]
            intsc = np.in1d(all_ions, SimIons)
            SimElmsIons = np.array(Ions)[intsc]
            SimElmsIons = [tuple(SimElmsIons[i]) for i in range(len(SimElmsIons))]
            spectra["SimIons"] = self.WrapSpectra(SimElmsIons, projectionSIW, sightparams, vel_mod=False, therm_mod=False)
            
            if self.ThermEff:
                spectra["SimIons_ThermalEffectsOff"] = self.WrapSpectra(SimElmsIons, projectionSIW, sightparams, vel_mod=False, therm_mod=True)

            if self.PecVelEff:
                spectra['SimIons_PeculiarVelocitiesOff'] = self.WrapSpectra(SimElmsIons, projectionSIW, sightparams, vel_mod=True, therm_mod=False)

            if self.ThermEff and self.PecVelEff:
                spectra['SimIons_ThermalPecVelOff'] = self.WrapSpectra(SimElmsIons, projectionSIW, sightparams, vel_mod=True, therm_mod=True)
        
        return spectra

    
    def WrapSpectra(self, Ions, projection, sightparams, vel_mod=False, therm_mod=False):
        ''' 
        Wraps and computes spectra for a given set of ions based on the projection and sightline parameters.

        Args:
            Ions (list): A list of tuples, each containing the element name and ion name (e.g., [('H', 'HI'), ('He', 'He II')]).
            projection (dict): The projection dictionary containing data for ion-weighted quantities such as densities,
                velocities, and temperatures for each ion.
            sightparams (dict): A dictionary containing the properties of the sightline such as velocity and pixel information.
            vel_mod (bool, optional): If True, modifies ion velocities to zero. Defaults to False.
            therm_mod (bool, optional): If True, modifies ion temperatures to a small value (0.1 K). Defaults to False.

        Returns:
            dict: A dictionary containing the spectra for each ion, including velocities, optical depths, densities,
                temperatures, and other ion information.

        Notes:
            This method computes the optical depth spectra for each ion based on the given projection and sightline parameters.
            It allows for modification of velocities and temperatures if specified.
        '''
        header = self.header
        spectra = {}
        vel_kms = sightparams['vel_kms']
        vunit = self.SetUnit(vardescription='Hubble velocity', Lunit=1e5, aFact=0, hFact=0)
        
        for ion in Ions:
            (element_name, ion_name) = ion
            weight = self.transitions[ion_name]["Mass"] * self.constants["amu"]
            lambda0 = self.transitions[ion_name]["lambda0"]
            f_value = self.transitions[ion_name]["f-value"]
            
            if lambda0 > 0:
                # Convert properties to CGS units
                nions = self.ToCGS(header, projection[ion_name]["Densities"]) / weight
                vions = self.ToCGS(header, projection[ion_name]["Velocities"]) / 1e5
                Tions = self.ToCGS(header, projection[ion_name]["Temperatures"])
                
                # Modify velocities and temperatures if desired
                if vel_mod:
                    vions = np.zeros_like(vions)
                if therm_mod:
                    Tions = np.zeros_like(Tions) + 0.1
                # Calculate the optical depth.
                spectrum = self.MakeOpticaldepth(
                    sightparams=sightparams,
                    weight=weight, 
                    lambda0=lambda0, 
                    f_value=f_value,
                    nions=nions, 
                    vions_kms=vions, 
                    Tions=Tions,
                    element_name=element_name
                )
                
                # Store computed spectra for each ion
                spectra[ion] = {
                    "Velocities": {'Value': vel_kms, 'Info': vunit},
                    "Optical depths": spectrum["Optical depths"],
                    "Densities": spectrum["Densities"],
                    "Velocities": spectrum["Velocities"],
                    "Temperatures": spectrum["Temperatures"],
                    "TotalIonColumnDensity": spectrum["TotalIonColumnDensity"],
                    "Mass": weight,
                    "lambda0": lambda0,
                    "f-value": f_value
                }
        
        return spectra

    
    def MakeOpticaldepth(self, sightparams={'sight_kms': 0.0, 'vel_kms': [0.0], 'pixel_kms': 1.0, 'pixel': 1.0},
                        weight=1.67382335232e-24, lambda0=1215.67, f_value=0.4164, 
                        nions=[0.0], vions_kms=[0.0], Tions=[0.0], element_name='Hydrogen'):
        ''' 
        Compute the optical depth for a given transition based on ionic density, temperature, and peculiar velocity.

        Args:
            sightparams (dict): A dictionary containing sightline parameters:
                - 'sight_kms' (float): Extent of the box in km/s.
                - 'vel_kms' (list of float): Array of velocity values in km/s.
                - 'pixel_kms' (float): Pixel size in km/s.
                - 'pixel' (float): Pixel size in cm.
            weight (float, optional): Mass of the ion in grams. Defaults to 1.67382335232e-24 (mass of hydrogen).
            lambda0 (float, optional): Rest wavelength of the transition in Angstroms. Defaults to 1215.67.
            f_value (float, optional): Oscillator strength of the transition. Defaults to 0.4164.
            nions (list of float, optional): Array of ionic densities in particles/cm^3. Defaults to [0.0].
            vions_kms (list of float, optional): Array of ionic velocities in km/s. Defaults to [0.0].
            Tions (list of float, optional): Array of ionic temperatures in K. Defaults to [0.0].
            element_name (str, optional): Name of the element (e.g., 'Hydrogen'). Defaults to 'Hydrogen'.

        Returns:
            dict: A dictionary containing:
                - 'Optical depths': Optical depth values and units.
                - 'Densities': Tau-weighted ion mass densities and units.
                - 'Velocities': Tau-weighted ion velocities and units.
                - 'Temperatures': Tau-weighted ion temperatures and units.
                - 'TotalIonColumnDensity': Total column density of ions and units.

        Notes:
            - The method computes the optical depth by integrating the column density of ions along the sightline.
            - It accounts for the b-parameter and includes the option for periodic boundary conditions.
            - Hubble velocity and velocity offsets are added to the ion velocities when computing the spectrum.
            - For hydrogen, the Lorentzian convolution is applied to the optical depth if Voigt profiles are enabled.
        '''
        box_kms = sightparams['sight_kms']
        vel_kms = sightparams['vel_kms']
        pixel_kms = sightparams['pixel_kms']
        pixel = sightparams['pixel']
        npix = len(vel_kms)
        tau = np.zeros_like(vel_kms)
        densities = np.zeros_like(vel_kms)
        velocities = np.zeros_like(vel_kms)  
        temperatures = np.zeros_like(vel_kms)

        # Set up the line properties and boundary conditions
        lines = Lines(v_kms=vel_kms, box_kms=box_kms, constants=self.constants, verbose=False, 
                    lambda0_AA=lambda0, f_value=f_value, periodic=self.periodic)
        
        # Convert from density to column density
        ioncolumns = nions * pixel  # in ions/cm^2
        total_column_density = np.sum(ioncolumns)  # in ions/cm^2
        dunit = self.SetUnit(vardescription="Total ion column density", Lunit=1.0, aFact=0.0, hFact=0.0)
        total_column_density = {'Value': total_column_density, "Info": dunit}

        # Compute b-parameter (velocity dispersion)
        bions_kms = np.sqrt(2 * self.constants["kB"] * Tions / weight) / 1e5
        
        # Add Hubble velocity and velocity offset to peculiar velocity
        vHubble_kms = box_kms * np.arange(len(vions_kms)) / len(vions_kms)
        voffset_kms = self.specparams['ODParams']['Veloffset']  # Default = 0 
        vions_tot_kms = vions_kms + vHubble_kms + voffset_kms
        
        spectrum = lines.gaussian(column_densities=ioncolumns, b_kms=bions_kms, vion_kms=vions_tot_kms, Tions=Tions)

        tau = spectrum['optical_depth']
        densities = spectrum['optical_depth_densities']
        velocities = spectrum['optical_depth_velocities']
        temperatures = spectrum['optical_depth_temperatures']
        pixel_velocities = spectrum['pixel_velocity_kms']

        # Apply Voigt for hydrogen profile convolution if required
        if not self.VoigtOff and element_name == "Hydrogen":
            tau = lines.convolvelorentz(tau)
        
        # Set units and convert mass density
        dunit = self.SetUnit(vardescription="Tau weighted ion mass density", Lunit=1.0, aFact=0.0, hFact=0.0)
        densities *= weight  # Convert to mass density
        densities = {'Value': densities, "Info": dunit}
        
        # Set units for velocities and temperatures
        vunit = self.SetUnit(vardescription="Tau weighted ion peculiar velocity", Lunit=1e5, aFact=0.0, hFact=0.0)
        velocities = {'Value': velocities, 'Info': vunit}
        tunit = self.SetUnit(vardescription="Tau weighted ion temperature", Lunit=1, aFact=0.0, hFact=0.0)
        temperatures = {'Value': temperatures, 'Info': tunit}
        tauunit = self.SetUnit(vardescription="Ionic optical depth", Lunit=1, aFact=0.0, hFact=0.0)
        tau = {'Value': tau, 'Info': tauunit}

        return {
            'Optical depths': tau, 
            'Densities': densities, 
            'Velocities': velocities, 
            'Temperatures': temperatures, 
            'TotalIonColumnDensity': total_column_density
        }

    def CGSunit(self, header, variable):
        ''' 
        Compute the conversion factor needed to convert simulation values to proper, h-free CGS units.

        Args:
            header (dict): Dictionary containing cosmological parameters of the simulation, 
                including 'Cosmo' with keys 'Redshift' and 'HubbleParam'.
            variable (dict): Dictionary containing information about the variable, 
                with 'Info' providing the scale exponents and conversion factor.

        Returns:
            float: The conversion factor to convert simulation values to CGS units.
                The formula is: CGSunit = CGSConversionFactor * h**hexpo * a**aexpo.
        '''
        # Dependence on expansion factor
        ascale = (1. / (1 + header["Cosmo"]["Redshift"]))**variable["Info"]["aexp-scale-exponent"]

        # Dependence on Hubble parameter
        hscale = (header["Cosmo"]["HubbleParam"])**variable["Info"]["h-scale-exponent"]
        
        return variable["Info"]["CGSConversionFactor"] * ascale * hscale


    def ToCGS(self, header, variable):
        ''' 
        Convert simulation values of the variable to proper CGS units (no h dependence).

        Args:
            header (dict): Dictionary containing cosmological parameters of the simulation.
            variable (dict): Dictionary containing the variable value and its information, 
                including the CGS conversion factor.

        Returns:
            float: The variable value converted to CGS units.
        '''
        return variable["Value"] * self.CGSunit(header, variable)


    def SetUnit(self, vardescription='text describing variable', Lunit=constants['Mpc'], aFact=1.0, hFact=1.0):
        ''' 
        Set the unit and conversion factors for a variable based on its length, scale, and Hubble parameter dependencies.

        Args:
            vardescription (str, optional): A description of the variable. Defaults to 'text describing variable'.
            Lunit (float, optional): The CGS conversion factor for the variable (e.g., length unit in cm). Defaults to constants['Mpc'].
            aFact (float, optional): The exponent factor for the expansion scale (a). Defaults to 1.0.
            hFact (float, optional): The exponent factor for the Hubble parameter (h). Defaults to 1.0.

        Returns:
            dict: A dictionary containing:
                - 'VarDescription': Description of the variable.
                - 'CGSConversionFactor': Conversion factor for the variable.
                - 'aexp-scale-exponent': Exponent for the scale factor a.
                - 'h-scale-exponent': Exponent for the Hubble parameter h.
        '''
        return {
            'VarDescription': vardescription, 
            'CGSConversionFactor': Lunit, 
            'aexp-scale-exponent': aFact, 
            'h-scale-exponent': hFact
        }
