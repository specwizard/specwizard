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


    def FileType(self, sim_type ='EAGLE', snap_type = 'LOS'):
        """
        This function will check if the simulation type and/or the simulation file is supported.
        If is the case will format this values into the dictionary with all the inputs. The input is not case-sensitive.  

        
        Args:
            sim_type (str): Simulation type, supported options ``'Eagle','Swift','Hydrangea'`` 
            snap_type (str): Type of output format full snapshots or line of sight. Valid options ``'Snapshot','LOS'``
                
        """

        snap = snap_type.lower()
        sym = sim_type.lower()
        sym_check = sym in self.sym_types
        snap_check = snap in self.snap_types

        if not sym_check:
            sys.exit("Simulation not supported")
        if not snap_check:
            sys.exit("snap_type must be a snapshot file or a line of sight (LOS) file.") 
            
        
        self.specparams['file_type']['sim_type']  = sym 
        self.specparams['file_type']['snap_type'] = snap 
    
        #return self.specparams 
    def SnapshotParams(self,path = '/cosma7/data/Eagle/ScienceRuns/Planck1/L0100N1504/PE/REFERENCE/data/los/', file = 'part_los_z3.275.hdf5'):
        """
        Will check if the path and file are valid. If so will be formatted into the specparams. 

        Args:
            path (str): /path/to/simulation/files
            file (str): simulationfile.hdf5
        """

        self.check_for_directory(path)

        self.check_hdf5(path+"/"+file)

        self.specparams['snapshot_params']['directory'] = path 
        self.specparams['snapshot_params']['file']      = file 
            
    def Sightline(self,nsight  = None, ProjectionAxes = None ,ProjectionStart = None, ProjectionLength=None,SightLength=None,ProjectionExtend=None,):

        """
        This functions generate the sightline characteristics. We recommend to see the example notebooks for further clarification.  If ``None``, the default values will be used.

        Args: 
            nsight (int): If using LOS simulation input, corresponds to the number of sightline in the file to read, is also used for labeling the sightline (``default=0``).
            ProjectionAxes Tuple[str]: Used for snapshot input. This create the relation between the projection axis and the simulation axis. The default case is when the projection axis is aligned with the simulation axis ``ProjectionAxes = ('projx','projy','projz') = ('simx','simy','simz')``. But if for example we want to do the projection along the $y-axis$ of the simulation then ``ProjectionAxes=('simx','simz','simy')``
            ProjectionStart Tuple[float]: The starting position as a fraction of the box-length of the sightline. If we take the default value, means that our sighline  will be centered in the simulation point $(x,z)=(0.5*\mathrm{boxsize},0.5*\mathrm{boxsize})$ and will start at the coordinate $z=0$. (``default=(0.5,0.5,0)``)
            ProjectionLength (float): Set as a fraction of the box size the length of the sightline. (``default=1``) 
            ProjectionExtend (float): Set as a fraction of the sightline length. This will insert the spectra into an array x times bigger than the sightline length. Is usefull for broadening effects. (``default=1.0``)
        """
        SimDict          = {'simx':0 , 'simy':1, 'simz':2}
        if ProjectionAxes == None:
            xaxis = yaxis = zaxis = None
        else:
            xaxis = SimDict[ProjectionAxes[0]]
            yaxis = SimDict[ProjectionAxes[1]]
            zaxis = SimDict[ProjectionAxes[2]]
        if ProjectionStart == None:
            xpos = ypos = zpos = None
        else:
            xpos,ypos,zpos = ProjectionStart
            
        sightline = {'ProjectionAxes':ProjectionAxes,
                        'x-axis': xaxis,
                        'y-axis': yaxis,
                        'z-axis': zaxis,
                        'ProjectionStart':ProjectionStart,
                        'x-position' :xpos,
                        'y-position' :ypos,
                        'z-position' :zpos,
                        'ProjectionLength': ProjectionLength,
                        'SightLength':SightLength,
                        'ProjectionExtend':ProjectionExtend,
                        'nsight':     nsight}
        
        if self.specparams['file_type']['snap_type'] == 'snapshot':
            sightline_default = {
                        'ProjectionAxes':('simx','simy','simz'),
                        'x-axis': 0,
                        'y-axis': 1,
                        'z-axis': 2,
                        'ProjectionStart':(0.5,0.5,0),
                        'x-position' :0.5,
                        'y-position' :0.5,
                        'z-position' :0,
                        'ProjectionLength': 1,
                        'SightLength':1,
                        'ProjectionExtend': {"extend":False,"extendfactor":3},
                        'nsight':     0}
                        
            default_keys = sightline_default.keys()

            for dkey in default_keys:
                if sightline[dkey]==None:
                    print("Warning! "+ dkey +" NOT found. Setting default value : "+str(sightline_default[dkey]))
                    sightline[dkey] = sightline_default[dkey]            
                    
        elif self.specparams['file_type']['snap_type'] == 'los': 
            if sightline['nsight']==None:
                print("Error! number of sightline (nsight) has not been assigned. \n This is necessary to select the LOS for the LOS file")
        
        self.specparams['sightline'] = sightline
        self.non_sense_checks(self.specparams)
        return self.specparams
        
    def SetIonTableParams(self,  table_type='specwizard_cloudy', iondir = '/cosma7/data/Eagle/SpecWizardCloudytables/HM12/', fname='', ions =  [('Hydrogen', 'H I') ],SFR_properties = {'modify_particle':True,'ignore_particle': True,'Temperature [K]':1e4},atomfile='./VpFit/atom.dat'):
        """
        Set the ion to use to calculate the spectra, the ionization tables to use for the calculation of the ion fractions and how to deal with star forming particles.

        Args: 
            table_type (str): Name of the ionizations tables to use current options ``'specwizard_cloudy'`` and ``'ploeckinger'``. The former corresponds to the Cloudy tables created from the Cloudy notebook based in HM2012, and the latter correspond to the tables described in Ploeckinger et al. (2020). (``default='specwizard_cloud'``)
            iondir (str): Directory path to the ionization table. (``default='HM12_table/path/in/cosma'``)
            fname (str): Name of ionization table, only needed for ploeckinger tables. (``default=''``)
            ions (list of tuples) : List of elements and Ions to do e,g ``[(Element1, Ion1, Ion2),(Element2, Ion1,...),...]``. ``default=[('Hydrogen', 'H I') ]``
            SFR_properties (dict) : if ``modify_particle: True`` we will deal with Star forming particles in one of two ways. If ``ignore particle = True``, will set to zero the IonizationFraction of SFR particles. If ``False`` we will set the Temperature of SFR by the indicated temperature value.  (default={'modify_particle':True,'ignore_particle': True,'Temperature [K]':1e4})  
        """
        ions_available = []
        table_types = ['specwizard_cloudy', 'ploeckinger']
        
        table_type = table_type.lower()
        if not table_type in table_types:
            print("Error: that is not a valid table type")
            sys.exit(-1)

        if table_type == 'specwizard_cloudy':
            try:
                iondir = iondir + '/'
                files  = os.listdir(iondir)
                ions_available   = []
                for file in files:
                    ions_available.append(file.split('.')[0])
            except:
                print("No files ionization tables found in directory {}".format(iondir))
                ions_available = []
        if table_type == 'Ploeckinger':
            try:
                file = iondir + '/' + fname
                files  = os.listdir(iondir)
                if not fname in files:
                    sys.exit(-1)
                hf       = h5py.File(file,  "r")
                elements = np.array(hf['ElementNamesShort'][...], dtype='str')
                hf.close()
            except:
                print("iondir = ", iondir)
                print("fname = ", fname)
                print("Doh!")
                sys.exit(-1)


        ion_properties = {'table_type':table_type, 
                            'iondir':iondir+'/', 
                            'fname':fname,
                            'ions-available':ions_available,
                            'Ions': ions,
                            'SFR_properties':SFR_properties,
                            'atomfile':atomfile}
        self.specparams['ionparams'] = ion_properties

        #set Ions and elements
        elements         = Elements(atomfile=atomfile)
        transitionparams = {}

        self.specparams["ionparams"]['IonizationBalance'] = IonTables(specparams  = self.specparams)

        #We always need hydrogen to be read since the ion fractions are calculated based on the hydrogen density 
        elements_check = [ions[i][0] for i in range(len(ions))]
        if "Hydrogen" not in elements_check:
            ions.insert(0,('Hydrogen','H I'))
        
        for element, ion in ions:
            pars = elements.ElementParameters(ElementNames=[element])
            transitionparams[ion] = {}
            transitionparams[ion]["Mass"]    = pars[element]["Weight"]

            try:
                transitionparams[ion]["lambda0"] = pars[element]["States"][ion]["Lines"]["Lambda0"][0]
                transitionparams[ion]["f-value"] = pars[element]["States"][ion]["Lines"]["f-value"][0]
                
            except:
                raise Exception("Unrecognized ion: "+ ion)
                
        self.specparams["ionparams"]["transitionparams"] = transitionparams    

        # Unique list of elemens
        ElementNames = []
        for (element, name) in ions:
            ElementNames.append(element)
        ElementNames    = np.unique(ElementNames)  # unique list of desired elements
        elementparams   = {}
        elementparams["ElementNames"] = ElementNames
        for element, ion in ions:
            pars = elements.ElementParameters(ElementNames=[element])
            elementparams[element] = {}
            elementparams[element]["Mass"] = pars[element]["Weight"]
        
        self.specparams["elementparams"]  = elementparams
        
    def SetODParams(self,VelOffset_kms = 0., PecVelEffectsOff = False, ThermalEffectsOff = False,VoigtOff = False):
        """
        User options that affect the calculation of the Optical Depth. 

        Args:
            ThermalEffectsOff (bool): Turn Off the thermal effecs in the calculation of OD's. (``default='False'``)
            PecVelEffectsOff (bool):  Turn off the effects of peculiar velocities. (``default='False'``)  
            VelOffset (float): Impose a shift in velocity space [km/s].  (``default=0.0``)
            VoigtOff (bool): If ``True``, will turn off the voigt profile for hydrogen. (``default='False'``)

        """
        
        
        self.specparams['ODParams']  = {'ThermalEffectsOff' : ThermalEffectsOff,
                                        'PecVelEffectsOff'  : PecVelEffectsOff, 
                                        'Veloffset'         : VelOffset_kms,
                                        'VoigtOff'          : VoigtOff }
     
    def SetLongSpectraParams(self,lambda_min = 945.,lambda_max = 8000., dlambda = 0.5, z_qsr = 3.0, delta_z = 0.01,file_dir = ''):
        """
        Set the parameters for the long spectra

        Args: 
            lambda_min (float): Set the minimum wavelength that will be included in the spectra. (``default = 945.``)
            lambda_max (float): Set the maximum wavelength that will be included in the spectra. (``default = 8000.``)
            dlambda (float): Set the pixel size of the spectrograph
            z_qsr (float): Redshift location of the "back light quasar". (``default = 3.0``)
            delta_z (float): Redshift tolerance. This will find the simulation file(s) that are within that tolerence for a given redshift section of the spectra. In case of finding multiplefiles one will be randomly selected. In case of not finding any file, the procedure will stop. (``default = 0.01``)
            file_dir (str):  Directory in which the simulation files that will constitute the longspectra are. SpecWizard will automatically read them and sort them by redshift. (``default=' '``) 
        """
        self.specparams['longspectra']  = {'lambda_min' : lambda_min,
                                        'lambda_max'    : lambda_max, 
                                        'dlambda'       : dlambda,
                                        'z_qsr'         : z_qsr ,
                                        'delta_z'       :delta_z,
                                        'file_dir'      :file_dir}
        
        
    def ExtraParams(self,periodic= True, kernel="Gauss",pixkms=1,atomfile='./VpFit/atom.dat',veloff=0,
                    ReadIonFrac = {'ReadIonFrac':False,
                        'ReadHydrogen':True,
                        'HI':"NeutralHydrogenAbundance",
                        'ReadHelium':False,
                        'He':"",
                        'fname_urchin':""} ):
        """
        This are parameters that are thought for specific or more advance users. 
        Args:
            periodic (bool) : Turn on/off periodic boundary condition effects (``default=True``).  
            Kernel  (string) : Selection of the SPH kernel for the interpolation in the line projection (``default=Bspline``). 
        """
        
        self.specparams['extra_parameters'] = { 'periodic': periodic,
                                                'Kernel': kernel ,
                                                'pixkms':  pixkms,
                                                'atomfile': atomfile,
                                                'Veloffset': veloff,
                                                'ReadIonFrac':ReadIonFrac}
        #self.non_sense_checks(self.specparams)

    def non_sense_checks(self,specparams):
        """
        This Fuction is set to check that the input makes sense.

        Args:
            specparams (dict)
        """
        sightline = specparams['sightline']
        extra_params = specparams['extra_parameters']
        #Checks if the axis provides is 0,1 and 2 (the sume must be 3)
        if (self.specparams['file_type']['snap_type'] == 'snapshot') and (sightline['x-axis']+sightline['y-axis']+sightline['z-axis'] != 3):
            raise Exception(" Invalid values for axis!")

        if (self.specparams['file_type']['snap_type'] == 'snapshot') and (extra_params['periodic']==True) and (sightline['ProjectionLength']<1):
            print("Warning! You are calculating a los shorter than the box-size \n and using periodic boundary conditions")

    def check_for_directory(self,path):
        """
        Checks if the path is a valid directory

        Args:
            path (str): Path to verify
        """
        bol_dir = os.path.isdir(path)
        if not bol_dir:
            raise Exception("Path: "+path +"\n is not a valid directory")
        return bol_dir

    def check_hdf5(self,fname):
        """
        Checks if the hdf5 file provided is a valid file.

        Args:
            fname (str): filename.hdf5
        """
        if not h5py.is_hdf5(fname):
            raise Exception(fname+" is NOT a valid hdf5 file")
            
    def read_from_yml(self,yml_file = 'Wizard.yml'):
        """
        This class read all the input data from a yml file 

        Args:
            yml_file (str): Name and path pointing to the parameter file. (default = "Wizard.yml")

        Returns:
            Formated dictionary with all the input that will be use for the rest of the program.
        """
        with open('Wizard.yml') as file:
            wizard_yml = yaml.load(file, Loader=yaml.FullLoader)
        
        # File type section from Wizard dictionary 
        sim_type   = wizard_yml['file_type']['sim_type']
        snap_type  = wizard_yml['file_type']['snap_type']

        self.FileType(sim_type=sim_type,snap_type=snap_type)

        # Snap type section from Wizard dictionary 
        snap_dir =  wizard_yml['snapshot_params']['directory']
        snap_file     =  wizard_yml['snapshot_params']['file'] 

        self.SnapshotParams(path=snap_dir,file=snap_file)

        # ion params section from Wizard dictionary 
        table_type = wizard_yml['ionparams']['table_type']
        iondir     = wizard_yml['ionparams']['iondir']
        fname      = wizard_yml['ionparams']['fname']
        SFR_properties = wizard_yml['ionparams']['SFR_properties']
        ions_array = np.array(wizard_yml['ionparams']['ions'])
        ions = [ (ions_array[i,0],ions_array[i,1]) for i in range(len(ions_array))]
        atomfile = wizard_yml['ionparams']['atomfile']

        self.SetIonTableParams(table_type=table_type, iondir = iondir,SFR_properties=SFR_properties , fname=fname, ions=ions,atomfile=atomfile)

        # ODParams  section from Wizard dictionary
        veloffset =  wizard_yml['ODParams']['VelOffset_kms']
        pecveloff =  wizard_yml['ODParams']['PecVelEffectsOff']
        thermalfoff = wizard_yml['ODParams']['ThermalEffectsOff']
        voigtoff    = wizard_yml['ODParams']['VoigtOff']
        self.SetODParams(VelOffset_kms=veloffset,PecVelEffectsOff=pecveloff ,ThermalEffectsOff=thermalfoff,VoigtOff=voigtoff)

        # longspectra 
        
        try:
            lambda_min  = wizard_yml['LongSpectra']['lambda_min']
            lambda_max  = wizard_yml['LongSpectra']['lambda_max']
            dlambda     = wizard_yml['LongSpectra']['dlambda']
            z_qsr       = wizard_yml['LongSpectra']['z_qsr']            
            delta_z     = wizard_yml['LongSpectra']['delta_z']
            file_dir    = wizard_yml['LongSpectra']['file_dir']

            self.SetLongSpectraParams(lambda_min = lambda_min,lambda_max = lambda_max, dlambda = dlambda, z_qsr = z_qsr, delta_z = delta_z,file_dir = file_dir)   
        except:
            pass 
            
        # extraparams
        try:
            
            periodic  = wizard_yml['extraparams']['periodic']
            pixelkms  = wizard_yml['extraparams']['pixkms']
            read_ion  = wizard_yml['extraparams']['ReadIonFrac']
            self.ExtraParams(periodic= True,pixkms=pixelkms,
                                ReadIonFrac = read_ion )
        
        except:
            pass
        #SightLine params 
        ProjectionAxes = wizard_yml['sightline']['ProjectionAxes']
        ProjectionStart   = wizard_yml['sightline']['ProjectionStart']
        projection_length = wizard_yml['sightline']['ProjectionLength']
        sightlength      = wizard_yml['sightline']['SightLength']
        projextend      = wizard_yml['sightline']['ProjectionExtend']
        nsight           = wizard_yml['sightline']['nsight']

        Wizard  = self.Sightline(nsight  = nsight, ProjectionAxes = (ProjectionAxes) ,ProjectionStart = (ProjectionStart), ProjectionLength=projection_length,SightLength=sightlength,ProjectionExtend=projextend)        
        return Wizard 
    def check_param(self,dic,param):
        try:
            dic[param]
            return True
        except:
            return False
