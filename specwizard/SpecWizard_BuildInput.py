# %load SpecWizard_BuildInput.py
#%%writefile SpecWizard_BuildInput.py 

import os
import h5py
import sys
import numpy as np
from SpecWizard_Elements import Elements
from SpecWizard_IonTables import IonTables
import yaml
import traceback


class Build_Input:
    """
    This class contains the relevant functions that assign and shape all the input that SpecWizard 
    need to work properly. This will be done by running the functions and making consistency checks. 
    If the checks are successful the input parameter will be formatted into the dictionary 'specparams'
    this dictionary will be use in the rest of the program.
    """
    def __init__(self):

        self.sym_types  = ['eagle','swift','hydrangea','colibre','illustries']
        self.snap_types = ["snapshot","los"]

        self.specparams = {}

        ReadIonFrac = {'ReadIonFrac':False,
                        'ReadHydrogen':True,
                        'HI':"NeutralHydrogenAbundance",
                        'ReadHelium':False,
                        'He':""}

        #we set extra parameters as default. This can be change by running the function self.extra_params().
        self.specparams['file_type']              = {}
        self.specparams['snapshot_params']        = {}
        self.specparams['ionparams']              = {}
        self.specparams['sightline']              = {}
        self.specparams['ODParams']  = {'ThermalEffectsOff' : False,
                                        'PecVelEffectsOff'  : False, 
                                        'Veloffset'         : 0,
                                        'VoigtOff'          : False}

        self.specparams['extra_parameters']       = {   'periodic': True,
                                                        'Kernel': "Bspline",
                                                        'pixkms': 1,
                                                        'Veloffset': 0,
                                                        'ReadIonFrac':ReadIonFrac}


    def FileType(self, sim_type='EAGLE', snap_type='LOS'):
        """
        This function checks if the simulation type and/or the simulation file is supported.
        If supported, it formats these values into the dictionary with all the inputs.
        The input is not case-sensitive.

        Args:
            sim_type (str): Simulation type, supported options 'Eagle', 'Swift', 'Hydrangea'
            snap_type (str): Type of output format, full snapshots or line of sight. Valid options 'Snapshot', 'LOS'
        """

        # Convert input to lowercase for case-insensitive comparison
        snap = snap_type.lower()
        sym = sim_type.lower()

        # Check if simulation type and snapshot type are supported
        sym_check = sym in self.sym_types
        snap_check = snap in self.snap_types

        # Exit if simulation type is not supported
        if not sym_check:
            sys.exit("Simulation not supported")

        # Exit if snapshot type is not valid
        if not snap_check:
            sys.exit("snap_type must be a snapshot file or a line of sight (LOS) file.")

        # Update dictionary with formatted values
        self.specparams['file_type']['sim_type'] = sym
        self.specparams['file_type']['snap_type'] = snap
    
    def SnapshotParams(self, path='/cosma7/data/Eagle/ScienceRuns/Planck1/L0100N1504/PE/REFERENCE/data/los/',
                       file='part_los_z3.275.hdf5'):
        """
        Checks if the path and file are valid. If so, they are formatted into specparams.

        Args:
            path (str): /path/to/simulation/files
            file (str): simulationfile.hdf5
        """

        # Check if the specified directory exists
        self.check_for_directory(path)

        # Check if the specified HDF5 file exists
        self.check_hdf5(path + "/" + file)

        # Update dictionary with formatted values
        self.specparams['snapshot_params']['directory'] = path
        self.specparams['snapshot_params']['file'] = file

            
    def Sightline(self, nsight=None, ProjectionAxes=None, ProjectionStart=None, ProjectionLength=None, SightLength=None,
                  ProjectionExtend=None):
        """
        Generates the sightline characteristics. We recommend seeing the example notebooks for further clarification.
        If None, the default values will be used.

        Args:
            nsight (int): If using LOS simulation input, corresponds to the number of sightlines in the file to read,
                          also used for labeling the sightline (default=0).
            ProjectionAxes (Tuple[str]): Used for snapshot input. This creates the relation between the projection axis
                                         and the simulation axis. The default case is when the projection axis is aligned
                                         with the simulation axis ('projx', 'projy', 'projz') = ('simx', 'simy', 'simz').
                                         But if, for example, we want to do the projection along the y-axis of the simulation,
                                         then ProjectionAxes=('simx', 'simz', 'simy').
            ProjectionStart (Tuple[float]): The starting position as a fraction of the box-length of the sightline.
                                           If we take the default value, it means that our sightline will be centered in the
                                           simulation point (x, z) = (0.5 * box size, 0.5 * box size) and will start at the coordinate z=0.
                                           (default=(0.5, 0.5, 0))
            ProjectionLength (float): Set as a fraction of the box size the length of the sightline. (default=1)
            ProjectionExtend (float): Set as a fraction of the sightline length. This will insert the spectra into an array x
                                      times bigger than the sightline length. Useful for broadening effects. (default=1.0)
        """

        # Dictionary for simulation axis mapping
        SimDict = {'simx': 0, 'simy': 1, 'simz': 2}

        # Check if ProjectionAxes and ProjectionStart are provided, if not, set to None
        if ProjectionAxes is None:
            xaxis = yaxis = zaxis = None
        else:
            # Map simulation axis based on ProjectionAxes
            xaxis = SimDict[ProjectionAxes[0]]
            yaxis = SimDict[ProjectionAxes[1]]
            zaxis = SimDict[ProjectionAxes[2]]

        # Check if ProjectionStart is provided, if not, set to None
        if ProjectionStart is None:
            xpos = ypos = zpos = None
        else:
            # Unpack ProjectionStart tuple
            xpos, ypos, zpos = ProjectionStart

        # Define sightline dictionary with input parameters
        sightline = {'ProjectionAxes': ProjectionAxes,
                     'x-axis': xaxis,
                     'y-axis': yaxis,
                     'z-axis': zaxis,
                     'ProjectionStart': ProjectionStart,
                     'x-position': xpos,
                     'y-position': ypos,
                     'z-position': zpos,
                     'ProjectionLength': ProjectionLength,
                     'SightLength': SightLength,
                     'ProjectionExtend': ProjectionExtend,
                     'nsight': nsight}

        # Check if the simulation type is 'snapshot'
        if self.specparams['file_type']['snap_type'] == 'snapshot':
            # Define default sightline parameters for snapshot
            sightline_default = {
                'ProjectionAxes': ('simx', 'simy', 'simz'),
                'x-axis': 0,
                'y-axis': 1,
                'z-axis': 2,
                'ProjectionStart': (0.5, 0.5, 0),
                'x-position': 0.5,
                'y-position': 0.5,
                'z-position': 0,
                'ProjectionLength': 1,
                'SightLength': 1,
                'ProjectionExtend': {"extend": False, "extendfactor": 3},
                'nsight': 0
            }

            # Get default keys
            default_keys = sightline_default.keys()

            # Set default values for missing parameters
            for dkey in default_keys:
                if sightline[dkey] is None:
                    print("Warning! " + dkey + " NOT found. Setting default value: " + str(sightline_default[dkey]))
                    sightline[dkey] = sightline_default[dkey]

        # Check if the simulation type is 'los'
        elif self.specparams['file_type']['snap_type'] == 'los':
            # Check if nsight is provided for LOS
            if sightline['nsight'] is None:
                print("Error! number of sightline (nsight) has not been assigned.\n"
                      "This is necessary to select the LOS for the LOS file")

        # Update specparams with sightline information
        self.specparams['sightline'] = sightline

        # Perform additional checks on the parameters
        self.non_sense_checks(self.specparams)

        # Return updated specparams
        return self.specparams

        
    def SetIonTableParams(self, table_type='specwizard_cloudy', iondir='/cosma7/data/Eagle/SpecWizardCloudytables/HM12/',
                          fname='', ions=[('Hydrogen', 'H I')],
                          SFR_properties={'modify_particle': True, 'ignore_particle': True, 'Temperature [K]': 1e4},
                          atomfile='./VpFit/atom.dat'):
        """
        Set the ion to use to calculate the spectra, the ionization tables to use for the calculation of the ion fractions,
        and how to deal with star forming particles.

        Args:
            table_type (str): Name of the ionizations tables to use, current options 'specwizard_cloudy' and 'ploeckinger'.
                             The former corresponds to the Cloudy tables created from the Cloudy notebook based on HM2012,
                             and the latter corresponds to the tables described in Ploeckinger et al. (2020).
                             (default='specwizard_cloud')
            iondir (str): Directory path to the ionization table. (default='HM12_table/path/in/cosma')
            fname (str): Name of ionization table, only needed for ploeckinger tables. (default='')
            ions (list of tuples): List of elements and ions to do, e.g., [('Element1', 'Ion1'), ('Element2', 'Ion1', ...), ...].
                                   (default=[('Hydrogen', 'H I')])
            SFR_properties (dict): If 'modify_particle: True', we will deal with star-forming particles in one of two ways.
                                   If 'ignore_particle = True', will set to zero the IonizationFraction of SFR particles.
                                   If False, we will set the Temperature of SFR by the indicated temperature value.
                                   (default={'modify_particle': True, 'ignore_particle': True, 'Temperature [K]': 1e4})
            atomfile (str): Path to the atom file. (default='./atom_dat.hdf5')
        """

        # Available ionization tables
        ions_available = []
        table_types = ['specwizard_cloudy', 'ploeckinger']

        # Convert table_type to lowercase for case-insensitive comparison
        table_type = table_type.lower()

        # Check if table_type is valid
        if not table_type in table_types:
            print("Error: that is not a valid table type")
            sys.exit(-1)

        # Check the availability of ionization tables based on the specified table_type
        if table_type == 'specwizard_cloudy':
            try:
                iondir = iondir + '/'
                files = os.listdir(iondir)
                ions_available = [file.split('.')[0] for file in files]
            except:
                print("No ionization tables found in directory {}".format(iondir))
                ions_available = []

        elif table_type == 'ploeckinger':
            try:
                file = iondir + '/' + fname
                files = os.listdir(iondir)
                if not fname in files:
                    sys.exit(-1)
                hf = h5py.File(file, "r")
                elements = np.array(hf['ElementNamesShort'][...], dtype='str')
                hf.close()
            except:
                print("iondir = ", iondir)
                print("fname = ", fname)
                print("Doh!")
                sys.exit(-1)

        # Set ion_properties dictionary
        ion_properties = {'table_type': table_type,
                          'iondir': iondir + '/',
                          'fname': fname,
                          'ions-available': ions_available,
                          'Ions': ions,
                          'SFR_properties': SFR_properties,
                          'atomfile': atomfile}
        self.specparams['ionparams'] = ion_properties

        # Set Ions and elements using Elements class
        elements = Elements(atomfile=atomfile)
        transitionparams = {}

        # Set IonizationBalance using IonTables class
        self.specparams["ionparams"]['IonizationBalance'] = IonTables(specparams=self.specparams)

        # Check if hydrogen is present in the specified ions, if not, insert it at the beginning
        elements_check = [ions[i][0] for i in range(len(ions))]
        if "Hydrogen" not in elements_check:
            ions.insert(0, ('Hydrogen', 'H I'))

        # Set transitionparams for each specified ion
        for element, ion in ions:
            pars = elements.ElementParameters(ElementNames=[element])
            transitionparams[ion] = {}
            transitionparams[ion]["Mass"] = pars[element]["Weight"]

            try:
                transitionparams[ion]["lambda0"] = pars[element]["States"][ion]["Lines"]["Lambda0"][0]
                transitionparams[ion]["f-value"] = pars[element]["States"][ion]["Lines"]["f-value"][0]
            except:
                raise Exception("Unrecognized ion: " + ion)

        # Update specparams with transitionparams
        self.specparams["ionparams"]["transitionparams"] = transitionparams

        # Unique list of elements
        ElementNames = np.unique([ion[0] for ion in ions])
        elementparams = {"ElementNames": ElementNames}

        # Set Mass for each specified element
        for element, ion in ions:
            pars = elements.ElementParameters(ElementNames=[element])
            elementparams[element] = {"Mass": pars[element]["Weight"]}

        # Update specparams with elementparams
        self.specparams["elementparams"] = elementparams


    def SetODParams(self, VelOffset_kms=0., PecVelEffectsOff=False, ThermalEffectsOff=False, VoigtOff=False):
        """
        User options that affect the calculation of the Optical Depth.

        Args:
            ThermalEffectsOff (bool): Turn Off the thermal effects in the calculation of OD's. (default='False')
            PecVelEffectsOff (bool):  Turn off the effects of peculiar velocities. (default='False')
            VelOffset_kms (float): Impose a shift in velocity space [km/s].  (default=0.0)
            VoigtOff (bool): If True, will turn off the Voigt profile for hydrogen. (default='False')
        """

        # Set ODParams dictionary with user options
        self.specparams['ODParams'] = {'ThermalEffectsOff': ThermalEffectsOff,
                                       'PecVelEffectsOff': PecVelEffectsOff,
                                       'Veloffset': VelOffset_kms,
                                       'VoigtOff': VoigtOff}

    def SetLongSpectraParams(self, lambda_min=945., lambda_max=8000., dlambda=0.5, z_qsr=3.0, delta_z=0.01, file_dir=''):
        """
        Set the parameters for the long spectra.

        Args:
            lambda_min (float): Set the minimum wavelength that will be included in the spectra. (default=945.0)
            lambda_max (float): Set the maximum wavelength that will be included in the spectra. (default=8000.0)
            dlambda (float): Set the pixel size of the spectrograph.
            z_qsr (float): Redshift location of the "backlight quasar". (default=3.0)
            delta_z (float): Redshift tolerance. This will find the simulation file(s) that are within that tolerance
                            for a given redshift section of the spectra. In case of finding multiple files,
                            one will be randomly selected. In case of not finding any file, the procedure will stop.
                            (default=0.01)
            file_dir (str): Directory in which the simulation files that will constitute the long spectra are.
                            SpecWizard will automatically read them and sort them by redshift. (default='')
        """
        # Set longspectra dictionary with user options
        self.specparams['longspectra'] = {'lambda_min': lambda_min,
                                          'lambda_max': lambda_max,
                                          'dlambda': dlambda,
                                          'z_qsr': z_qsr,
                                          'delta_z': delta_z,
                                          'file_dir': file_dir}

    def ExtraParams(self, periodic=True, kernel="Gauss", pixkms=1, veloff=0,
                    ReadIonFrac={'ReadIonFrac': False,
                                  'ReadHydrogen': True,
                                  'HI': "NeutralHydrogenAbundance",
                                  'ReadHelium': False,
                                  'He': "",
                                  'fname_urchin': ""}):
        """
        Parameters for specific or more advanced users.

        Args:
            periodic (bool): Turn on/off periodic boundary condition effects (default=True).
            Kernel (string): Selection of the SPH kernel for the interpolation in the line projection (default=Bspline).
            pixkms (int): Pixel size in km/s.
            veloff (float): Velocity offset. (default=0)
            ReadIonFrac (dict): Dictionary for reading ion fractions. (default={'ReadIonFrac': False,
                                                                                  'ReadHydrogen': True,
                                                                                  'HI': "NeutralHydrogenAbundance",
                                                                                  'ReadHelium': False,
                                                                                  'He': "",
                                                                                  'fname_urchin': ""})
        """
        # Set extra_parameters dictionary with user options
        self.specparams['extra_parameters'] = {'periodic': periodic,
                                               'Kernel': kernel,
                                               'pixkms': pixkms,
                                               'Veloffset': veloff,
                                               'ReadIonFrac': ReadIonFrac}


    def non_sense_checks(self, specparams):
        """
        Check if the input makes sense.

        Args:
            specparams (dict): Dictionary of spectrograph parameters.
        """
        sightline = specparams['sightline']
        extra_params = specparams['extra_parameters']

        # Check if the axis provides is 0, 1, and 2 (the sum must be 3)
        if (self.specparams['file_type']['snap_type'] == 'snapshot') and (sightline['x-axis'] + sightline['y-axis'] + sightline['z-axis'] != 3):
            raise Exception("Invalid values for axis!")

        # Check if using periodic boundary conditions and calculating a LOS shorter than the box-size
        if (self.specparams['file_type']['snap_type'] == 'snapshot') and (extra_params['periodic'] is True) and (
                sightline['ProjectionLength'] < 1):
            print("Warning! You are calculating a LOS shorter than the box-size \n and using periodic boundary conditions")

    def check_for_directory(self, path):
        """
        Check if the path is a valid directory.

        Args:
            path (str): Path to verify.
        """
        bol_dir = os.path.isdir(path)
        if not bol_dir:
            raise Exception("Path: " + path + "\n is not a valid directory")
        return bol_dir

    def check_hdf5(self, fname):
        """
        Check if the HDF5 file provided is a valid file.

        Args:
            fname (str): Filename.hdf5
        """
        if not h5py.is_hdf5(fname):
            raise Exception(fname + " is NOT a valid HDF5 file")
            
    def read_from_yml(self, yml_file='Wizard.yml'):
        """
        Read all the input data from a YAML file.

        Args:
            yml_file (str): Name and path pointing to the parameter file. (default="Wizard.yml")

        Returns:
            Formatted dictionary with all the input that will be used for the rest of the program.
        """
        with open(yml_file) as file:
            wizard_yml = yaml.load(file, Loader=yaml.FullLoader)

        # File type section from Wizard dictionary
        sim_type = wizard_yml['file_type']['sim_type']
        snap_type = wizard_yml['file_type']['snap_type']
        self.FileType(sim_type=sim_type, snap_type=snap_type)

        # SnapshotParams section from Wizard dictionary
        snap_dir = wizard_yml['snapshot_params']['directory']
        snap_file = wizard_yml['snapshot_params']['file']
        self.SnapshotParams(path=snap_dir, file=snap_file)

        # IonParams section from Wizard dictionary
        table_type = wizard_yml['ionparams']['table_type']
        iondir = wizard_yml['ionparams']['iondir']
        fname = wizard_yml['ionparams']['fname']
        SFR_properties = wizard_yml['ionparams']['SFR_properties']
        ions_array = np.array(wizard_yml['ionparams']['ions'])
        ions = [(ions_array[i, 0], ions_array[i, 1]) for i in range(len(ions_array))]
        atomfile = wizard_yml['ionparams']['atomfile']
        self.SetIonTableParams(table_type=table_type, iondir=iondir, SFR_properties=SFR_properties, fname=fname,
                               ions=ions, atomfile=atomfile)

        # ODParams section from Wizard dictionary
        veloffset = wizard_yml['ODParams']['VelOffset_kms']
        pecveloff = wizard_yml['ODParams']['PecVelEffectsOff']
        thermalfoff = wizard_yml['ODParams']['ThermalEffectsOff']
        voigtoff = wizard_yml['ODParams']['VoigtOff']
        self.SetODParams(VelOffset_kms=veloffset, PecVelEffectsOff=pecveloff, ThermalEffectsOff=thermalfoff,
                         VoigtOff=voigtoff)

        # LongSpectra section
        try:
            lambda_min = wizard_yml['LongSpectra']['lambda_min']
            lambda_max = wizard_yml['LongSpectra']['lambda_max']
            dlambda = wizard_yml['LongSpectra']['dlambda']
            z_qsr = wizard_yml['LongSpectra']['z_qsr']
            delta_z = wizard_yml['LongSpectra']['delta_z']
            file_dir = wizard_yml['LongSpectra']['file_dir']
            self.SetLongSpectraParams(lambda_min=lambda_min, lambda_max=lambda_max, dlambda=dlambda, z_qsr=z_qsr,
                                      delta_z=delta_z, file_dir=file_dir)
        except:
            pass

        # ExtraParams section
        try:
            periodic = wizard_yml['extraparams']['periodic']
            pixelkms = wizard_yml['extraparams']['pixkms']
            read_ion = wizard_yml['extraparams']['ReadIonFrac']
            self.ExtraParams(periodic=True, pixkms=pixelkms, ReadIonFrac=read_ion)
        except:
            pass

        # SightLine params
        ProjectionAxes = wizard_yml['sightline']['ProjectionAxes']
        ProjectionStart = wizard_yml['sightline']['ProjectionStart']
        projection_length = wizard_yml['sightline']['ProjectionLength']
        sightlength = wizard_yml['sightline']['SightLength']
        projextend = wizard_yml['sightline']['ProjectionExtend']
        nsight = wizard_yml['sightline']['nsight']
        Wizard = self.Sightline(nsight=nsight, ProjectionAxes=ProjectionAxes, ProjectionStart=ProjectionStart,
                                ProjectionLength=projection_length, SightLength=sightlength, ProjectionExtend=projextend)
        return Wizard

    def check_param(self, dic, param):
        """
        Check if a parameter exists in a dictionary.

        Args:
            dic (dict): Dictionary to check.
            param (str): Parameter to check.

        Returns:
            bool: True if the parameter exists, False otherwise.
        """
        try:
            dic[param]
            return True
        except:
            return False
