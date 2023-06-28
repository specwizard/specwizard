
import os
import h5py
import sys
import numpy as np
from SpecWizard_Elements import Elements
from SpecWizard_IonTables_test import IonTables


class Build_Input:
    '''
    This class contains the relevant functions that assign and shape all the input that SpecWizard 
    need to work properly. This will be done by running the functions and making consistency checks. 
    If the checks are sucesesfull the inputparameter will be formated into the dictionary 'specparams'
    this dictionary will be use in the rest of the program.
    '''
    
    def __init__(self):
    
        self.sym_types  = ['eagle','swift','hydrangea','colibre','illustries']
        self.snap_types = ["snapshot","los"]
        
        self.specparams = {}
        
        ReadIonFrac = {'ReadIonFrac':True,
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
        '''
        This function will check if the simulation type and/or the simulation file is supported. \n If this is the case will format this values into specparams. The format is not case-sensitive.  
        Input: sim_type: String, simulation type, supported options ['Eagle','Swift','Hydrangea'] 
               snap_type: String, type of output format full snapshots or line of sight  figths, valid options ['Snapshot','LOS']
               
       '''

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
    
#        return self.specparams 
    def SnapshotParams(self,path = '/cosma7/data/Eagle/ScienceRuns/Planck1/L0100N1504/PE/REFERENCE/data/los/', file = 'part_los_z3.275.hdf5'):
        
        '''
        Will check if the path and file are valid. If so will be formated into the specparams.
        Inputs:
        path: String, /path/to/simulation/files
        file: String, simulationfile.hdf5
        '''
        
        self.check_for_directory(path)
        
        self.check_hdf5(path+"/"+file)
        
        self.specparams['snapshot_params']['directory'] = path 
        self.specparams['snapshot_params']['file']      = file 
        
    ProjectionAxes   = (Projx,Projy,Projz) = ('simy','simz','simx')
    ProjectionStart  = (0.1 , 0.5, 0 )
    ProjectionLength = 1
    def Sightline(self,nsight  = None, ProjectionAxes = None ,ProjectionStart = None, ProjectionLength=None,SightLength=None,ProjectionExtend=None,):

        '''
        This generate the sightline characteristics:
        xwiz: Int, default=0, set which coordinate of the simulation data will be set as the x-axis.
        ywiz: Int, default=1, set which coordinate of the simulation data will be set as the y-axis.
        zwiz: Int, default=2, set which coordinate of the simulation data will be set as the z-axis.
        xpos : Float, default=0.5, set as a fraction of the box size the x-position of the sightline.
        ypos : Float, default=0.5, set as a fraction of the box size the y-position of the sightline.
        zpos : Float, default=0, set as a fraction of the box size where in the z coordinate the sightline starts.
        los-light: Float, default=1, set as a fraction of the box size the lenght of the sightline.
        nsight: Int, default=0, If using LOS simulation, corresponds to the number of sightline in the file to read, is used for labeling the sightline.

        '''
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
        '''
        Set the ion to use to calculate the spectra, the ionization tables to use for the calculation of the ion fractions and how to deal with star forming particles.
        Input:
        table_type: String, default='specwizard_cloud', name of the ionizations tables to use current options 'specwizard_cloudy' and 'ploeckinger'. The former corresponds to the Cloudy tables created from the Cloudy notebook based in HM2012, and the latter correspond to the tables described in Ploeckinger et al. (2020).
        iondir: String, default='HM12 table path in cosma', path for the ionization table.
        fname: String, default='',name of ionization tabel, only needed for ploeckinger tables. 
        ions: List of tuples, default=[('Hydrogen', 'H I') ], List of elements and Ions to do e,g [(Element1, Ion1, Ion2),(Element2, Ion1,...),...]
        SFR_properties, Dictionary, default={'modify_particle':True,'ignore_particle': True,'Temperature [K]':1e4}, if modiy_particle: True we will deal with Star forming particles in one of two ways. If ignore particle = True, will set to zero the IonizationFraction of SFR particles. If False we will set the Temperature of SFR by the indicated temperature value.   
        '''
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
                print("Warning: Default values for transition parameters are in use.")
                transitionparams[ion]["lambda0"] = 1000
                transitionparams[ion]["f-value"] = 1e-5
                
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
        '''
            User options that affect the calculation of the Optical Dept. 
            Input:
            ThermalEffectsOff: Bolean, default='False', Turn Off the thermal effecs in the calculation of OD's.
            PecVelEffectsOff: Bolean, default='False',  Turn off the effects of peculiar velocities .
            VelOffset: Float, default='0.0',Impose a shift in . 
            VoigtOff: Bolean, default='False', If True, will turn off the voigt profile for hydrogen

        '''
        
        
        self.specparams['ODParams']  = {'ThermalEffectsOff' : ThermalEffectsOff,
                                        'PecVelEffectsOff'  : PecVelEffectsOff, 
                                        'Veloffset'         : VelOffset_kms,
                                        'VoigtOff'          : VoigtOff}
        
        
    def ExtraParams(self,periodic= True, kernel="Gauss",pixkms=1,atomfile='./VpFit/atom.dat',veloff=0,
                    ReadIonFrac = {'ReadIonFrac':True,
                       'ReadHydrogen':True,
                       'HI':"NeutralHydrogenAbundance",
                       'ReadHelium':False,
                       'He':""} ):
        '''
        This are parameters that are thought for specific or more advance users. 
        Inputs:
            periodic: Bolean, default=True. Turn on/off periodic boundary condition effects. 
            Kernel  : string, default=Bspline. Selection of the SPH kernel for the interpolation in the line projection. 
            Veloff  : Float[kms],  default=0. Will add this value to all the peculiar velocities. Units km/s. 
        '''
        
        
            
        
        self.specparams['extra_parameters'] = { 'periodic': periodic,
                                                'Kernel': kernel ,
                                                'pixkms':  pixkms,
                                                'atomfile': atomfile,
                                                'Veloffset': veloff,
                                                'ReadIonFrac':ReadIonFrac}
       # self.non_sense_checks(self.specparams)

    def non_sense_checks(self,specparams):
        sightline = specparams['sightline']
        extra_params = specparams['extra_parameters']
        
        if (self.specparams['file_type']['snap_type'] == 'snapshot') and (sightline['x-axis']+sightline['y-axis']+sightline['z-axis'] != 3):
            print("ERROR! Invalid values for axis!")

        if (self.specparams['file_type']['snap_type'] == 'snapshot') and (extra_params['periodic']==True) and (sightline['ProjectionLength']<1):
            print("Warning! You are calculating a los shorter than the box-size \n and using periodic boundary conditions")

    def check_for_directory(self,path):
        bol_dir = os.path.isdir(path)
        if not bol_dir:
            sys.exit("Path: "+path +"\n is not a valid directory")
        return bol_dir
    
    def check_hdf5(self,fname):
        if not h5py.is_hdf5(fname):
            sys.exit("Not a valid hdf5 file")
            
            
    def check_param(self,dic,param):
        try:
            dic[param]
            return True
        except:
            return False