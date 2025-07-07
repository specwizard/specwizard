import h5py 
import glob
from .Phys import ReadPhys
import re
import unyt
import numpy as np
import swiftsimio as sw
import pyread_eagle as read_eagle
import hydrangea as hy
from .SimulationInputKeys import get_simkeys
import copy
from .SpecWizard_Elements import Elements

# physical constants in cgs units
constants = ReadPhys()

class InputFunctions:
    
    def __init__(self,fileparams={},header=None):
        """"
        This class containes a set of functions that are commonly shared among the different classes that read the data from simulations.
        
        Parameters
        ----------
        
        fileparams: dictionary
        The Wizard dictionary created by SpecWizard_BuildInput
        """ 
        
    
        self.fileparams = fileparams
        self.fdir       = fileparams["snapshot_params"]["directory"]
        self.fname      = self.fdir + '/' + fileparams["snapshot_params"]["file"]
        self.simtype    = fileparams["file_type"]["sim_type"]
        self.snaptype   = fileparams["file_type"]["snap_type"]
        self.readIonFrac = fileparams['extra_parameters']['ReadIonFrac']['ReadIonFrac'] 
        sim_keys  = get_simkeys(self.simtype)
        groupdic  = sim_keys[self.snaptype]
        self.groupdic = groupdic
        self.sightline = self.fileparams['sightline']
        groupname = groupdic['groupname']
        self.groupname = groupname
        
        if header != None:
            self.header = header
        
    def read_group(self, groupname ='Header',alt_fname=None):
        """
        Read the attributes of a specific hdf5 data group
        
        Parameters
        ----------
        
        groupname : str
        Name of the data group to read attributes.
        
        
        Returns
        -------
        
        out : dictionary 
        Dictionary with the name (key) and value of the attribute. 
        """
    # read all entries for this particular hdf5 group
        if alt_fname == None:
            fname = self.fname
        else:
            fname = alt_fname
        hfile = h5py.File(fname, "r")
        group    = hfile[groupname]
        grp_dict = {}
        for k in group.attrs.keys():
            grp_dict[k]= group.attrs[k]

        hfile.close()

        return dict(grp_dict)
    
    def set_unit(self, vardescription = 'text describing variable', Lunit=constants['Mpc'], aFact=1.0, hFact=1.0):
        """
        Formats into a dictionary the information for converting to CGS units. Using the toCGS function.  
        
        Parameters
        ----------
        
        vardescription : str
        Description of the variable e.g Mass of the gas particle
        
        Lunit  : float 
        Conversion factor to CGS.
        
        aFact : float
        Exponent of the scale factor $a$ for conversion to physical units. 
        
        hFact : float
        Exponent of the Hubble parameter $h=H/100$, some quantities in some simulations are divided by h. 
        
        Returns
        -------
        
        out : dictionary
        Formated dictionary with the necesary informations to conver to physical CGS units. 
        """
        return {'VarDescription': vardescription, 'CGSConversionFactor':Lunit, 'aexp-scale-exponent' :aFact, 'h-scale-exponent': hFact}

    
    def Hubble(self):
        """
         Hubble constant at redshift z, in units 1/s
        """
        Hz  = self.header["Cosmo"]["H0"]
        z   = self.header["Cosmo"]["Redshift"]
        Hz *= np.sqrt(self.header["Cosmo"]["OmegaMatter"] * (1.+z)**3 + self.header["Cosmo"]["OmegaLambda"])
        return Hz
    
    def FormatTxt(self,text):
        """
        Format text into lower case and eliminate spaces
        
        Parameters
        ----------
        
        text : str
        
        Returns
        -------
        
        out : str
        
        Text formated with lower case and no white spaces
        """
        return text.lower().replace(" ","")
    
    def set_sightlineparams(self,sightline,header):
        """
        Adds to the Wizard['sightline'] dictionary extra information about of the sightline, This is use specifically if a LOS file is used.
        
        Parameters
        ----------
        
        sightline : dictionary 
        sightline dictionary from the Wizard dictionary
        
        header    : dictionary 
        Header read from the simulation data. 
        
        Returns
        -------
        
        out : dictionary
        Returns the updated sightline dictionary/
        """
        
        self.header = header
        if self.snaptype == 'los':
            # We read extra information about the sightline located in the data file. Such position of the sightline, orientation, boxsize ...
            groupdic     = self.groupdic
            groupname    = groupdic['groupname'].format(sightline['nsight'])
            sightline_s  = self.read_group(groupname = groupname)
            try:
                num_particles = sightline_s[groupdic['Number_of_part_this_los']]
            except: 
                num_particles = 0
            box_size = header["BoxSize"]['Value'][0]
            xaxis  = sightline_s[groupdic['x-axis']]
            yaxis  = sightline_s[groupdic['y-axis']]
            zaxis  = sightline_s[groupdic['z-axis']]
            xpos   = sightline_s[groupdic['x-position']] / box_size.value
            ypos   = sightline_s[groupdic['y-position']] / box_size.value
            zpos   = 0 

            sightline    = { 'nsight': sightline['nsight'], 
                            'Number_of_part_this_los': num_particles,
                             'x-axis': xaxis,
                             'x-position': xpos,
                             'y-axis': yaxis,
                             'y-position': ypos,
                             'z-axis': zaxis,
                             'z-position': zpos}
            if (self.simtype == 'swift') or (self.simtype == 'colibre'):
                sightline['z-position'] = 0.
                sightline['x-axis']     = xaxis[0]
                sightline['y-axis']     = yaxis[0]
                sightline['z-axis']     = zaxis[0]  
            sightline['short-LOS'] = False
            sightline['ProjectionLength'] = 1 
            sightline['ProjectionExtend'] = {"extend":False,"extendfactor":1}

        # compute the velocity extent of the simulation volume in the sightline direction
        z_axis  = sightline["z-axis"]        
        box_cgs = self.to_physical(header["BoxSize"])[z_axis].in_cgs()       # in proper cm
        box_kms = box_cgs.to('km') * self.Hubble()# in km/s
        unit    = self.set_unit(vardescription='Extent of simulation volume in direction of projection in terms of Hubble velocity',
                        Lunit  = 1,
                        aFact  = 0.,
                        hFact  = 0.)
        
        sightline["Boxkms"] = {'Value': box_kms, 'Info': unit}
        sightline["Boxkms"]['Info'].pop('CGSConversionFactor','None')
        #
        box      = header["BoxSize"]["Value"][z_axis]
        box      = {'Value':box, 'Info': header["BoxSize"]['Info']}
        sightline["Box"] = box
        
        sightline['short-LOS'] = False
        if sightline['ProjectionLength'] < 1:
            sightline['short-LOS'] = True

        return sightline
        
    def BuildLOS(self,sightline,header):
        """
        Build LOS from full simulation box. Takes the spatial information provided by the user during BuildInput, and shape it so it can be used by     ReadEagle, swiftsim io. This to build the desire sightline.  
        
        Parameters
        ---------- 
        
        sightline : dictionary
        sightline dictionary from Wizard['sightline'] 
        
        header    : dictionary 
        Header read from the simulation data. 
        
        Returns
        -------
        (tuple): tuple containing:
        
            axis_dic (dictionary) Maps the user specified sightline axes to the simulation axes. 
            xpos (float) x-position of the sightline in the same units of length from the simulation
            ypos (float) y-position of the sightline in the same units of length from the simulation
            zpos (float) z-position of the sightline in the same units of length from the simulation
            los_length-(float) sightline length in the same units of length from the simulation.
            msl_x3 (float) three times the mean separation length between the particles. 
        
        """
        
        def MeanSeparationLenght(N,V):
            
            """
            We calculate the mean separation length to define the region that will be masked by the simulation reading function.

            Parameters
            ----------

            V : float
                Volume of the simulation
            N : int
                Total number of particles

            Returns
            -------

            out : float 
                Mean separation length. 
            """
            ninv = V/N
            return pow(ninv,1/3)
                    
        # We calculate the mean separation length to define the region that will be masked by the simulation reading function. 
        N = header['NumPartTot']['Value']
        V = np.prod(header['BoxSize']['Value'])
        msl_x3 = 3*MeanSeparationLenght(N,V)
        # We get the los_lenght as a fraction and multiply by the box size so we can extract particles from 
        # z-position to z-position+los_length 
        
        if sightline['ProjectionLength'] < 1:
            sightline['short-LOS'] = True

        los_length  = sightline['ProjectionLength']        
        xpos        = sightline['x-position']
        ypos        = sightline['y-position']
        zpos        = sightline['z-position']
        
        box         = header['BoxSize']['Value'][0]
    
        los_length  *= box
        xpos        *= box
        ypos        *= box
        zpos        *= box
        
        xmin,xmax = [xpos-msl_x3,xpos+msl_x3]
        ymin,ymax = [ypos-msl_x3,ypos+msl_x3]
        zmin,zmax = [zpos,zpos+los_length]
        axis_dic = {sightline['x-axis']:[xmin,xmax],
                    sightline['y-axis']:[ymin,ymax],
                    sightline['z-axis']:[zmin,zmax]}

        #Since the simulation reading functions use the simulation convention of axis we adjust them.
        sim_xmin,sim_xmax = axis_dic[0]
        sim_ymin,sim_ymax = axis_dic[1]
        sim_zmin,sim_zmax = axis_dic[2]
        
        
        return axis_dic,xpos,ypos,zpos,los_length,msl_x3
    
    
    
    def ReadAndShapeIonFrac(self,read_variable,particles,groupname):
        """
        
         Reads and formats the ion fraction information from simulations that have that information. Currently only supported by Colibre
         
         Parameters
         ----------
         
         read_variable : function
         Function that reads the particle data for the simulation. This is function is defined in SpecWizard_ReadData
         
         particles : dictionary
         Particle dictionary containg all the other particle properties. This is for proper conversion of the ion fractions.
         
         groupname : str
         Name of the hdf5 group in which the ion fractions are located. 
         
         Returns
         -------
         
         out : dictionary
         Formated dictionary with the ion fractions.
         
        """
        
        groupdic         = self.groupdic
        userions         = self.fileparams['ionparams']['Ions']

#        if self.simtype == 'eagle':
            #In this case we use the output from urchin the ion fractions are in different files from the snapshot but they can be read using readeagle
#                IonFracs = {}

#                IonFracs[ion] = {'Value': IonFrac['Value'],'Info': IonFrac['Info']}

        if self.simtype == 'colibre':

            # In the case of colibre the ion fraction is saved as a fraction of hydrogen content, while specwizard uses ion as a fraction of the element of the ion. 
            atomfile = self.fileparams['ionparams']['atomfile']
            elements = Elements(atomfile)
            elementnames = [name[0] for name in userions] 
            elements_info = elements.ElementParameters(ElementNames=elementnames)

            IonFracs = {}
            colfile    = h5py.File(self.fname,'r')
            col_ions   = colfile['SubgridScheme']['NamedColumns']['SpeciesFractions'][...].astype('str')

#            print((1./(1+self.header["Cosmo"]["Redshift"]))**particles['Densities']['Info']["aexp-scale-exponent"])
            handy_dens = copy.deepcopy(particles['Densities'])
            handy_dens['Value'] = self.assing_unit_unyt(particles['Densities'],'Densities')
            Abundances = particles['Abundances']
            dens_cgs   = self.to_physical(handy_dens).in_cgs()
            nH         = dens_cgs * Abundances['Hydrogen']['Value'] / constants['mH']

            col_ions_formated = np.array([self.FormatTxt(col_ion) for col_ion in col_ions])
            req_ions_formated = np.array([self.FormatTxt(ion) for _, ion in userions])
            ions_we_have = np.where(np.intersect1d(col_ions_formated,req_ions_formated))[0]
            
            for elements, ion in np.array(userions)[ions_we_have]:
                mass_e   =  elements_info[elements]['Weight'] * constants['amu']
                ne       = dens_cgs * Abundances[elements]['Value'] / mass_e
                ion_formated = self.FormatTxt(ion)
                ions_indx = np.where(col_ions_formated==ion_formated)[0]
                ionname   = str(col_ions[ions_indx][0])
                IonFrac =  read_variable(varname = groupname + '/' + groupdic['IonFractions'] + '/'+ ionname)
                IonFrac['Value'] = IonFrac['Value'] * (nH/ne)
                
                IonFracs[ion] = {'Value': IonFrac['Value'],'Info': IonFrac['Info']}
        return IonFracs    
    
    def set_fractions_metallicity(self,read_variable,particles):
        """
        Reads the element fraction and metallicities of particles. The metallicy is used in some ionization tables (e.g Ploekinger)
        
        Since simulations do not contain metallicity or element fraction. This fuction attempts to read the simulation data, if it fails it set some default values
        Sets metallicities to zero 
        Sets elemen fraction of hydrogen and helium to primordial values.
        
        Parameters
        ----------
        
        read_variable : function
        Function that reads the particle data for the simulation. This is function is defined in SpecWizard_ReadData

        particles : dictionary
        Particle dictionary containg all the other particle properties. 
        
        
        Returns
        -------
        
        out : tuple containing:
            element2do (list) : list of strings of the intersection between the user input elements and the elements available in the ionization table. 
            abundances (dictionary) : formated dictionary with the available element abundances.
            metallicity (dictionary) : formated dictionary with the particle metallicities. 
            
        
        
        """
        groupdic         = self.groupdic 
        groupname        = self.groupname
        elementnames     = groupdic['elementnames']
        userions         = self.fileparams['ionparams']['Ions']
        elements_we_want = np.array([userions[i][0] for i in range(len(userions))])
        elements2do      = np.intersect1d(elements_we_want,elementnames)
        abundances = {}
        unit       = self.set_unit(vardescription = 'mass fraction of element',Lunit=1.0, aFact=0.0, hFact=0.0) 

        if self.snaptype == 'los':
            groupname    = groupdic['groupname'].format(self.sightline['nsight'])

        
        try:    
        

            if len(elements2do) == 0:
                print("ERROR! No valid element or element name!")


            for elementname in elements2do:

                values  = read_variable(varname = groupname +'/'+ groupdic['ElementAbundance']+'/'+ elementname)['Value']
                info    = unit
                abundances[elementname] = {'Value':values, 'Info': unit}
        except:
            print("Element fraction not found using primordial quantities...")
            hydrogen    = np.zeros_like(particles['Densities']['Value']) + 1. - constants['Y']
            abundances['Hydrogen'] = {'Value':hydrogen, 'Info':unit}
            helium      = np.zeros_like(particles['Densities']['Value']) + constants['Y']
            abundances['Helium']  = {'Value':helium, 'Info':unit}

        #Check if metallicities are available. 
        try: 
            metallicity  =read_variable(varname = groupname + '/' + groupdic['Metallicities'])
        except:
            print("Warning! Metallicities not found. Setting them to zero.")
            metallicity = np.zeros_like(particles['Densities']['Value'])
            unit        = self.set_unit(vardescription = 'mass fraction of metals', 
                                       Lunit=1.0, aFact=0, hFact=0)
            metallicity = {'Value':metallicity, 'Info':unit} 
            

        
        return elements2do,abundances,metallicity

    
    def flip_los(self,data,sightinfo):
        """
        This function will put the long projected axis as the third column in the output. So the longest part of the sightline will always be at [:,2]
        
        Parameters
        ----------
        
        data : array 
        Data array contaning 3D data (e,g positions, velocities)
        
        sightinfo : dictionary
        Dictionary with the sightline info this should be part of the Wizard dictionary from BuildInput
        
        Returns
        -------
        
        out : array
        Formated data array. With the modified axis. 
        
        """
        
        temp = data.copy()
        temp[:,0] = data[:,sightinfo['x-axis']]
        temp[:,1] = data[:,sightinfo['y-axis']]
        temp[:,2] = data[:,sightinfo['z-axis']]

        return temp
    
    
    def ToCGS(self, variable):
        """ 
        Return simulations values for this variable in proper cgs units (no h)
        
        Parameters
        ----------
        
        variable : dictionary 
        Dictionary of attribute to be converted to CGS.
        
        Returns
        -------
        
        out : array 
        array of the CGS converted quantity. 
        
        """
        return variable["Value"] * self.CGSunit(variable)
    
    
    
    def CGSunit(self, variable):
        """ 
        Use the information in the variable to compute the factor needed to convert
        simulation values to proper, h-free cgs units.
        This is of the form:
        proper value = simulation value * CGSunit, where
        CGSunit = CGSConversionFactor * h**hexpo * a**aexpo
        This is used by the ToCGS function. 
        
        Parameters
        ----------
        
        variable : dictonary
        Formated dictionary of the attribute we wish to covert to CGS.
        
        Returns
        ------- 
        
        out : float
        Conversion factor for physical CGS. 
        
        """
        #header     = self.read_group()
        #dependence on expansion factor
        #ascale     = (1./(1+self.header["Cosmo"]["Redshift"]))**variable["Info"]["aexp-scale-exponent"]

        # dependence on hubble parameter
        #hscale     = (self.header["Cosmo"]["HubbleParam"])**variable["Info"]["h-scale-exponent"]
        
        #
        return variable["Info"]["CGSConversionFactor"] #* ascale * hscale
    

    def set_SFR(self,read_variable,groupname,particles):
        """
        Attempts to read star formation rates, if not creates a zero array. For some users dealing with star forming particles in important and how to deal with them is part of BuildInput. 
        
        Parameters
        ----------
        
        read_variable : function
        Function that reads the particle data for the simulation. This is function is defined in SpecWizard_ReadData

        particles : dictionary
        Particle dictionary containg all the other particle properties. 

        groupname : str
        Name of the hdf5 group in which the SFR is located. 
        
        
        """
        groupdic     = self.groupdic 
        particlesSFR = {}
        try:

            SFR    = read_variable(varname = groupname + '/' + groupdic['StarFormationRate'])
            particlesSFR['StarFormationRate'] = SFR
        except:
            print('Warning! Not able to read Star formation Rate properties...setting SFR to zero!')
            SFR         = np.zeros_like(particles['Masses']['Value']) 
            unit        = self.set_unit(vardescription = 'Star Formation Rate', 
                           Lunit=1.0, aFact=0, hFact=0)
            particlesSFR['StarFormationRate'] = {'Value':SFR, 'Info':unit} 
            
        return particlesSFR['StarFormationRate']
    def unyt_dict(self,):
        '''
        returns unyt dictionary 
        '''
        unyt_dict = {}

        unyt_dict["In"] = {"Positions": unyt.cm,
                    "Masses":unyt.g,
                    "Velocities":unyt.cm/unyt.s,
                    "Densities":unyt.g/unyt.cm**3,
                    "Temperatures":unyt.K,
                    "SmoothingLengths":unyt.cm,
                    "Abundances":unyt.dimensionless,
                    "SimulationIonFractions":unyt.dimensionless,
                    "Metallicities":unyt.dimensionless,
                    "StarFormationRate":unyt.g/unyt.s,
                    "BoxSize":unyt.cm}


        unyt_dict["Out"] = {"Positions": unyt.Mpc,
                    "Masses":1e10*unyt.Msun,
                    "Velocities":unyt.km/unyt.s,
                    "Densities":1e10*unyt.Msun/unyt.Mpc**3,
                    "Temperatures":unyt.K,
                    "SmoothingLengths":unyt.Mpc,
                    "Abundances":unyt.dimensionless,
                    "SimulationIonFractions":unyt.dimensionless,
                    "Metallicities":unyt.dimensionless,
                    "StarFormationRate":unyt.Msun/unyt.s,
                    "BoxSize":unyt.Mpc}

        return unyt_dict 

    def assing_unit_unyt(self,quantity,item_name):

        '''
        Converts the quantity into a unyt array. 
        '''

        unyt_dictionary = self.unyt_dict()
        item_value = unyt.unyt_array(self.ToCGS(quantity),unyt_dictionary["In"][item_name])
        return item_value.to(unyt_dictionary["Out"][item_name])



    def set_particle_data_to_unyt(self,part_data):
        
        test_unyt = {}


        for item in part_data.keys():
            if item == "Elements":
                continue
            if item == "Abundances" or item == "SimulationIonFractions":
                test_unyt[item] = {}

                for abun in part_data[item]:

                    test_unyt[item][abun] = {}
                    test_unyt[item][abun]['Value'] = self.assing_unit_unyt(part_data[item][abun],item)
                    test_unyt[item][abun]['Info']  =  {}
                    test_unyt[item][abun]['Info']['VarDescription'] = item+" from the simulation particle data. The units are in co-moving and units are descrived in the unyt array."             
                    test_unyt[item][abun]['Info']['aexp-scale-exponent'] = part_data[item][abun]['Info']['aexp-scale-exponent']
                    test_unyt[item][abun]['Info']['h-scale-exponent'] = part_data[item][abun]['Info']['h-scale-exponent']
                continue

            else:
                test_unyt[item] = {}
                test_unyt[item]['Value'] = self.assing_unit_unyt(part_data[item],item)
                test_unyt[item]['Info']  =  {}
                test_unyt[item]['Info']['VarDescription'] = item+" from the simulation particle data. The units are in co-moving and units are descrived in the unyt array."             
                test_unyt[item]['Info']['aexp-scale-exponent'] = part_data[item]['Info']['aexp-scale-exponent']
                test_unyt[item]['Info']['h-scale-exponent'] = part_data[item]['Info']['h-scale-exponent']

        return test_unyt

    def to_physical(self,variable):
        '''
        Convert unyt arrya from comoving to physical
        '''

        ascale     = (1./(1+self.header["Cosmo"]["Redshift"]))**variable["Info"]["aexp-scale-exponent"]

        # dependence on hubble parameter
        hscale     = (self.header["Cosmo"]["HubbleParam"])**variable["Info"]["h-scale-exponent"]


        return variable['Value']*ascale*hscale

        


    
    
    
class ReadEagle:

    
    def __init__(self,fileparams={}):
        
        """
        Class that contains functions to read particle data from the eagle simulation. It uses the python port of ReadEagle:
        https://github.com/kyleaoman/pyread_eagle
        
        Parameters
        ----------
        
        
        fileparams : dictionary
        The Wizard dictionary, output from SpecWizard_BuildInput
        """
        self.fileparams = fileparams
        self.fdir       = fileparams["snapshot_params"]["directory"]
        self.fname      = self.fdir + '/' + fileparams["snapshot_params"]["file"]
        self.simtype    = fileparams["file_type"]["sim_type"]
        self.snaptype   = fileparams["file_type"]["snap_type"]
        self.readIonFrac = fileparams['extra_parameters']['ReadIonFrac']['ReadIonFrac'] 
        sim_keys        = get_simkeys(self.simtype)
        groupdic        =  sim_keys[self.snaptype]
        self.groupdic = groupdic
        
        groupname = groupdic['groupname']
        self.groupname = groupname        
        self.inputfunc = InputFunctions(fileparams)
    
    
    
    def read_header(self):
        """
        Reads header from the eagle simulation. 
        
        
        """

        header = self.inputfunc.read_group()
 
        # the unit and h-dependence of Eagle is not stated; We assume it is in cMpc/h
        boxsize  = np.array([1.0, 1.0, 1.0 ]) * header['BoxSize']
        boxunit  = self.inputfunc.set_unit(vardescription="Extent of simulation volume", 
                                Lunit=constants['Mpc'], 
                                aFact=1.0, 
                                hFact=-1.0)
        box      = {'Value':boxsize, 'Info':boxunit}
        boxsize  = self.inputfunc.assing_unit_unyt(box,'BoxSize')
        box      = {'Value':boxsize, 'Info':boxunit}
        box['Info'].pop('CGSConversionFactor',None)

        
        h        = header['HubbleParam']
        #
        cosmo    = {'Redshift'    : header['Redshift'], 
                    'HubbleParam' : h,
                    'OmegaMatter' : header['Omega0'],
                    'OmegaBaryon' : header['OmegaBaryon'],
                    'OmegaLambda' : header['OmegaLambda']}

        numpartval   = header['NumPart_Total'][0]
        numpartunit  = self.inputfunc.set_unit(vardescription="Total number of gas particles in the simulation", 
                                Lunit=1, 
                                aFact=0, 
                                hFact=0)
        numpart      = {'Value':numpartval,'Info':numpartunit}
        numpart['Info'].pop('CGSConversionFactor',None)

        initMassTable = header['MassTable']

        DMmass        = initMassTable[1]
        OmegaDM       = (cosmo['OmegaMatter']-cosmo['OmegaBaryon'])
        Gasmass       = DMmass * (cosmo['OmegaBaryon']/OmegaDM)

        Dmunit        = self.inputfunc.set_unit(vardescription="Initial Dark Matter mass", 
                                Lunit=1.989e43,                                     
                                aFact=0, 
                                hFact=-1)
        Gmunit        = self.inputfunc.set_unit(vardescription="Initial Gass particle mass", 
                                Lunit=1.989e43,                                     
                                aFact=0, 
                                hFact=-1)
        DMmass        = {'Value':DMmass,'Info':Dmunit}
        DMmass  = self.inputfunc.assing_unit_unyt(DMmass,'Masses')
        DMmass        = {'Value':DMmass,'Info':Dmunit}
        DMmass['Info'].pop('CGSConversionFactor',None)

        Gasmass       = {'Value':Gasmass,'Info':Gmunit}
        Gasmass  = self.inputfunc.assing_unit_unyt(Gasmass,'Masses')
        Gasmass       = {'Value':Gasmass,'Info':Gmunit}
        Gasmass['Info'].pop('CGSConversionFactor',None)
        MassTable     = {'DarkMatterMass':DMmass,'GasMass':Gasmass}
        
        # compute some extra variables
        H0            = h * 100 * 1e5 * unyt.cm / unyt.s / constants["Mpc"]      # H0 in 1/s
        rhoc          = 3*H0**2 / (8*np.pi*constants["Grav"]) # critical density in g/cm^3
        cosmo["H0"]   = H0
        cosmo["rhoc"] = rhoc
        cosmo["rhob"] = rhoc * cosmo["OmegaBaryon"]

        Header = {'BoxSize':box, 'Cosmo' : cosmo, 'NumPartTot':numpart, 'MassTable':MassTable}

        return Header

    def read_particles(self):
        """
        Read the all the needed particle properties for the calculation of the spectra. Both for snapshots and line of sight files. 
        """
        
        groupname = self.groupname
        groupdic  = self.groupdic
        header    = self.read_header() 
        sightline = self.fileparams['sightline']
        
        sightline = self.inputfunc.set_sightlineparams(sightline,header)
        
        if self.snaptype == 'los':
            groupname    = groupdic['groupname'].format(sightline['nsight'])

        if self.snaptype == 'snapshot':
            
            
            sightline
            axis_dic,xpos,ypos,zpos,los_length,msl_x3  = self.inputfunc.BuildLOS(sightline,header)
            sim_xmin,sim_xmax = axis_dic[0]
            sim_ymin,sim_ymax = axis_dic[1]
            sim_zmin,sim_zmax = axis_dic[2]
        
        
            self.RE_snap = read_eagle.EagleSnapshot(self.fname)
            self.RE_snap.select_region(sim_xmin,sim_xmax,sim_ymin,sim_ymax,sim_zmin,sim_zmax)            
            SmoothingL = self.RE_snap.read_dataset(0,'SmoothingLength') * unyt.Mpc
            Positions  = self.RE_snap.read_dataset(0,'Coordinates') * unyt.Mpc


            mask_x     =  ((Positions[:,sightline['x-axis']]>(xpos-SmoothingL)) & (Positions[:,sightline['x-axis']]<(xpos+SmoothingL)))
            mask_y     =  ((Positions[:,sightline['y-axis']]>(ypos-SmoothingL)) & (Positions[:,sightline['y-axis']]<((ypos+SmoothingL))))
            mask_z     =  ((Positions[:,sightline['z-axis']]>zpos) & (Positions[:,sightline['z-axis']]<((zpos+los_length))))

            self.los_mask = np.where((mask_x) & (mask_y) & (mask_z))[0]

            
            
        particles = {'Densities'        : self.read_variable(varname = groupname + '/' + groupdic['Densities']),
                     'SmoothingLengths' : self.read_variable(varname = groupname + '/' + groupdic['SmoothingLengths']),
                     'Masses'           : self.read_variable(varname = groupname + '/' + groupdic['Masses']),
                     'Positions'        : self.read_variable(varname = groupname + '/' + groupdic['Positions']),
                     'Velocities'       : self.read_variable(varname = groupname + '/' + groupdic['Velocities']),
                     'Temperatures'     : self.read_variable(varname = groupname + '/' + groupdic['Temperatures']),
                     'Elements'         : None, 
                     'Abundances'       : None  }
            
        elements,abundances,metallicities    = self.inputfunc.set_fractions_metallicity(self.read_variable,particles)    
        
        particles['Elements']   = elements
        particles['Abundances'] = abundances
        particles['Metallicities'] = metallicities
        
        if self.fileparams['ionparams']['SFR_properties']['modify_particle']:            
            particles['StarFormationRate']   = self.inputfunc.set_SFR(self.read_variable,groupname,particles)

        if self.readIonFrac:
            IonFracs = {}
            fname_urchin = self.fileparams['extra_parameters']['ReadIonFrac']['fname_urchin']
            h1_name    = self.fileparams['extra_parameters']['ReadIonFrac']['HI']
            self.inputfunc = InputFunctions(self.fileparams)
            self.RE_snap = read_eagle.EagleSnapshot(fname_urchin)
            self.RE_snap.select_region(sim_xmin,sim_xmax,sim_ymin,sim_ymax,sim_zmin,sim_zmax)
            IonFrac = self.read_variable(varname = groupname + '/' + h1_name,alt_fname=fname_urchin)        
            IonFracs["H I"] = {'Value': IonFrac['Value'],'Info': IonFrac['Info']}
            particles['SimulationIonFractions'] = IonFracs

        particles['Positions']['Value'] = self.inputfunc.flip_los(particles['Positions']['Value'],sightline)
        particles['Velocities']['Value'] = self.inputfunc.flip_los(particles['Velocities']['Value'],sightline)        
        
        return particles,sightline
        
        
    def read_variable(self,varname='PartType0/Density',alt_fname=None):
        """
        Reads a particular hdf5 group and formats it's attributes
        
        Parameters
        ----------
        
        varname : str
        name of the hdf5 group with full hdf5 path
        """
        
        if self.snaptype == 'los':
            info   = self.inputfunc.read_group(groupname = varname)
            # The meta data of the velocities is incorrect in gadget los files
            # The output in the file is the internal velocity unit, a^2 dx/dt.
            # To correct this to proper peculiar velocity, we need to divide by the expansion factor
            if self.snaptype == 'los':
                if self.IsVelocity(varname):
                    info["VarDescription"]      = 'Co-moving velocities. Physical v_p = a dx/dt  = Velocities a^(-1) U_V [cm/s]'
                    info['aexp-scale-exponent'] = -1

                hfile  = h5py.File(self.fname, "r")
                values = np.array(hfile[varname][...])
                hfile.close()
        
        if self.snaptype =='snapshot':
            varname_frmtd = varname.replace("PartType0/","")
            values = self.RE_snap.read_dataset(0,varname_frmtd)[self.los_mask]
            if alt_fname == None:
                info   = self.inputfunc.read_group(groupname = varname)
            else:
                info   = self.inputfunc.read_group(groupname = varname,alt_fname=alt_fname)
        return {'Value': values, 'Info': info}
        
    def IsVelocity(self, name):
        """
        Check whether we are reading the velocity variable. This check is needed since line of sight files from eagle uses a wrong conversion factors for velocities. 
        """
        words    = name.split('/')
        Velocity = False
        for word in words:
                if (word == 'Velocity') or (word == 'Velocities'):
                    Velocity = True
        return Velocity


    
class ReadSwift:
    
    """
    Class that contains functions to read particle data produced by the SWIFT code. it uses the python library swiftsim.io 
    https://github.com/SWIFTSIM/swiftsimio

    Parameters
    ----------


    fileparams : dictionary
        The Wizard dictionary, output from SpecWizard_BuildInput
    """
    
    
    def __init__(self,fileparams={}):
        self.fileparams = fileparams
        self.fdir       = fileparams["snapshot_params"]["directory"]
        self.fname      = self.fdir + '/' + fileparams["snapshot_params"]["file"]
        self.simtype    = fileparams["file_type"]["sim_type"]
        self.snaptype   = fileparams["file_type"]["snap_type"]
        self.readIonFrac = fileparams['extra_parameters']['ReadIonFrac']['ReadIonFrac'] 
        sim_keys        = get_simkeys(self.simtype)
        groupdic        =  sim_keys[self.snaptype]
        self.groupdic = groupdic
        
        groupname = groupdic['groupname']
        self.groupname = groupname        
        self.inputfunc = InputFunctions(fileparams)
    
    
    
    def read_header(self):
        
        """
        Reads header as produced by SWIFT. This includes the COLIBRE simulation.  
        
        
        """
        
        header    = self.inputfunc.read_group()
        cosmology = self.inputfunc.read_group('Cosmology')
        units     = self.inputfunc.read_group('Units')
        boxsize   = np.array(header['BoxSize'])
        boxunit   = self.inputfunc.set_unit(vardescription="Extent of simulation volume", Lunit=units['Unit length in cgs (U_L)'][0], aFact=1.0, hFact=0.0) # co-moving, no-h
        box       = {'Value':boxsize, 'Info':boxunit}
        boxsize  = self.inputfunc.assing_unit_unyt(box,'BoxSize')
        box      = {'Value':boxsize, 'Info':boxunit}
        box['Info'].pop('CGSConversionFactor',None)
        
        #
        H0        = cosmology['H0 [internal units]'][0] / units['Unit time in cgs (U_t)'][0] # in 1/s
        H0        = H0 * constants['Mpc'] / 1e5     # in km/s/Mpc
        h         = H0 / 100                        # in units 100km/s/Mpc

        cosmo     = {'Redshift'   : cosmology['Redshift'][0], 
                    'HubbleParam' : h,
                    'OmegaMatter' : cosmology['Omega_m'][0],
                    'OmegaBaryon' : cosmology['Omega_b'][0],
                    'OmegaLambda' : cosmology['Omega_lambda'][0]}

        numpartval   = header['NumPart_Total'][0]
        numpartunit  = self.inputfunc.set_unit(vardescription="Total number of gas particles in the simulation", 
                                Lunit=1, 
                                aFact=0, 
                                hFact=0)
        numpart      = {'Value':numpartval,'Info':numpartunit}  
        numpart['Info'].pop('CGSConversionFactor',None)


        initMassTable = header['InitialMassTable']

        DMmass        = initMassTable[1]
        Gasmass       = initMassTable[0]

        Dmunit        = self.inputfunc.set_unit(vardescription="Initial Dark Matter mass", 
                                Lunit=1.98841e43,                                     
                                aFact=0, 
                                hFact=0)
        Gmunit        = self.inputfunc.set_unit(vardescription="Initial Gass particle mass", 
                                Lunit=1.98841e43,                                     
                                aFact=0, 
                                hFact=0)
        DMmass        = {'Value':DMmass,'Info':Dmunit}
        DMmass  = self.inputfunc.assing_unit_unyt(DMmass,'Masses')
        DMmass        = {'Value':DMmass,'Info':Dmunit}

        Gasmass       = {'Value':Gasmass,'Info':Gmunit}

        Gasmass  = self.inputfunc.assing_unit_unyt(Gasmass,'Masses')
        Gasmass       = {'Value':Gasmass,'Info':Gmunit}
        Gasmass['Info'].pop('CGSConversionFactor',None)

        MassTable     = {'DarkMatterMass':DMmass,'GasMass':Gasmass}
        
        
        # compute some extra variables
        H0            = h * 100 * 1e5 * unyt.cm / unyt.s / constants["Mpc"]      # H0 in 1/s
        rhoc          = 3*H0**2 / (8*np.pi*constants["Grav"]) # critical density in g/cm^3
        cosmo["H0"]   = H0
        cosmo["rhoc"] = rhoc
        cosmo["rhob"] = rhoc * cosmo["OmegaBaryon"]
            
        Header = {'BoxSize':box, 'Cosmo' : cosmo, 'NumPartTot':numpart, 'MassTable':MassTable}
        
        return Header
    
    
    
    def read_particles(self):
        """
        Read the all the needed particle properties for the calculation of the spectra. Both for snapshots and line of sight files. 
        """        
        
        groupname = self.groupname
        groupdic  = self.groupdic
        header    = self.read_header() 
        sightline = self.fileparams['sightline']
        
        sightline = self.inputfunc.set_sightlineparams(sightline,header)
        
        if self.snaptype == 'los':
            groupname    = groupdic['groupname'].format(sightline['nsight'])

        if self.snaptype == 'snapshot':
            
            
            sightline
            axis_dic,xpos,ypos,zpos,los_length,msl_x3  = self.inputfunc.BuildLOS(sightline,header)
            sim_xmin,sim_xmax = axis_dic[0]
            sim_ymin,sim_ymax = axis_dic[1]
            sim_zmin,sim_zmax = axis_dic[2]
        
        
            swif_mask = sw.mask(self.fname)
            Sboxsize = swif_mask.metadata.boxsize[0]
            oneMpc =  Sboxsize/Sboxsize.value                   #Since swiftsimIO uses units we have to 
            load_region = [[sim_xmin.value*oneMpc,sim_xmax.value*oneMpc],[sim_ymin.value*oneMpc,sim_ymax.value*oneMpc],[sim_zmin.value*oneMpc,sim_zmax.value*oneMpc]]
            swif_mask.constrain_spatial(load_region)

            self.SW_snap = sw.load(self.fname, mask=swif_mask)
            SmoothingL = self.SW_snap.gas.smoothing_lengths
            Positions  = self.SW_snap.gas.coordinates


            mask_x     =  ((Positions[:,sightline['x-axis']]>(-SmoothingL+xpos)) & (Positions[:,sightline['x-axis']]<(SmoothingL+xpos)))
            mask_y     =  ((Positions[:,sightline['y-axis']]>(-SmoothingL+ypos)) & (Positions[:,sightline['y-axis']]<((SmoothingL+ypos))))
            mask_z     =  ((Positions[:,sightline['z-axis']]>zpos) & (Positions[:,sightline['z-axis']]<((los_length+zpos))))

            self.los_mask = np.where((mask_x) & (mask_y) & (mask_z))[0]

            
            
        particles = {'Densities'        : self.read_variable(varname = groupname + '/' + groupdic['Densities']),
                     'SmoothingLengths' : self.read_variable(varname = groupname + '/' + groupdic['SmoothingLengths']),
                     'Masses'           : self.read_variable(varname = groupname + '/' + groupdic['Masses']),
                     'Positions'        : self.read_variable(varname = groupname + '/' + groupdic['Positions']),
                     'Velocities'       : self.read_variable(varname = groupname + '/' + groupdic['Velocities']),
                     'Temperatures'     : self.read_variable(varname = groupname + '/' + groupdic['Temperatures']),
                     'Elements'         : None, 
                     'Abundances'       : None  }
            
        elements,abundances,metallicity    = self.inputfunc.set_fractions_metallicity(self.read_variable,particles)    
        
        particles['Elements']   = elements
        particles['Abundances'] = abundances
        particles['Metallicities'] = metallicity
        
        
        FWHM = 0.362
        particles['SmoothingLengths']["Value"] /= FWHM
        print("We divide Swift's smoothing length by {0:1.3f} to convert from FWHM to extent of finite support".format(FWHM))

        if (self.simtype == 'swift' and self.readIonFrac):
            print("this is happening")
            field_name= self.fileparams['extra_parameters']['ReadIonFrac']['HI']
            ionfrac = self.read_variable(varname = groupname + '/' + field_name)
            particles['SimulationIonFractions'] = {}
            particles['SimulationIonFractions']["H I"] = ionfrac
        if (self.simtype == 'colibre'):
            
            particles['SimulationIonFractions'] = self.inputfunc.ReadAndShapeIonFrac(self.read_variable,particles,groupname)
        
        if self.fileparams['ionparams']['SFR_properties']['modify_particle']:
            
            particles['StarFormationRate']   = self.inputfunc.set_SFR(self.read_variable,groupname,particles)

        particles['Positions']['Value'] = self.inputfunc.flip_los(particles['Positions']['Value'],sightline)
        particles['Velocities']['Value'] = self.inputfunc.flip_los(particles['Velocities']['Value'],sightline)        
        
        return particles,sightline
    
    
    def read_variable(self,varname='PartType0/Density'):
        """
        Reads a particular hdf5 group and formats it's attributes
        
        Parameters
        ----------
        
        varname : str
        name of the hdf5 group with full hdf5 path
        """
        
        if self.snaptype == 'los':
            hfile  = h5py.File(self.fname, "r")
            values = np.array(hfile[varname][...])
            hfile.close()
            info_s  = self.inputfunc.read_group(groupname = varname)

        if self.snaptype =='snapshot':
            
            varname_frmtd = varname.split("/")
            
            swiftsimio_format = self.swiftswimio_format(varname = varname_frmtd[1])

            if varname_frmtd[1] == self.groupdic['ElementAbundance']:
                element_format = self.swiftswimio_format(varname = varname_frmtd[2])
                values = (getattr(getattr(self.SW_snap.gas,swiftsimio_format),element_format).value)[self.los_mask]
                info_s   = self.inputfunc.read_group(groupname = varname_frmtd[0]+'/'+varname_frmtd[1])

            elif varname_frmtd[1] == self.groupdic['IonFractions']:
                ion_format = varname_frmtd[2]
                values = (getattr(getattr(self.SW_snap.gas,swiftsimio_format),ion_format).value)[self.los_mask]
                info_s   = self.inputfunc.read_group(groupname = varname_frmtd[0]+'/'+varname_frmtd[1])            
            else:
                
                values = (getattr(self.SW_snap.gas,swiftsimio_format).value)[self.los_mask]
                info_s   = self.inputfunc.read_group(groupname = varname)
                
        info    = {'VarDescription': info_s['Description'] + info_s['Expression for physical CGS units'],
                           'CGSConversionFactor': info_s['Conversion factor to CGS (not including cosmological corrections)'][0],
                           'h-scale-exponent': info_s['h-scale exponent'][0],
                           'aexp-scale-exponent': info_s['a-scale exponent'][0]}
        return {'Value': values, 'Info': info}
    
  

    def swiftswimio_format(self,varname="Coordinates"):
        
        """
        Formats the string so it can be used for swiftsimio 
        
        Parameters
        ----------
        
        varname : str
        Name of the particle property to read from swift using swiftsim.io
        
        Returns
        -------
        
        out : str
        Formated string 
        """
        varname_formated = re.sub(r"([A-Z])", r"_\1", varname)

        return varname_formated[1:].lower()

    
class ReadHydrangea:
    """
    Class that contains functions to read particle data produced by the hydrangea simulation. it uses the python library to read hydrangea
    https://github.com/SWIFTSIM/swiftsimio

    """
    
    def __init__(self,fileparams={}):
        self.fileparams = fileparams
        self.fdir       = fileparams["snapshot_params"]["directory"]
        self.fname      = self.fdir + '/' + fileparams["snapshot_params"]["file"]
        self.simtype    = fileparams["file_type"]["sim_type"]
        self.snaptype   = fileparams["file_type"]["snap_type"]
        self.readIonFrac = fileparams['extra_parameters']['ReadIonFrac']['ReadIonFrac'] 
        sim_keys        = get_simkeys(self.simtype)
        groupdic        =  sim_keys[self.snaptype]
        self.groupdic = groupdic
        
        groupname = groupdic['groupname']
        self.groupname = groupname        
        self.inputfunc = InputFunctions(fileparams)
    
    
    
    def read_header(self):

             # read header information and store in default format
        header = self.inputfunc.read_group()

        # the unit and h-dependence of Eagle is not stated; We assume it is in cMpc/h
        boxsize  = np.array([1., 1.0, 1.0 ]) * header['BoxSize']
        boxunit  = self.inputfunc.set_unit(vardescription="Extent of simulation volume", 
                                Lunit=constants['Mpc'], 
                                aFact=1.0, 
                                hFact=-1.0)
        box      = {'Value':boxsize, 'Info':boxunit}
        boxsize  = self.inputfunc.assing_unit_unyt(box,'BoxSize')
        box      = {'Value':boxsize, 'Info':boxunit}
        box['Info'].pop('CGSConversionFactor',None)

        h        = header['HubbleParam']
        #
        cosmo    = {'Redshift'    : header['Redshift'], 
                    'HubbleParam' : h,
                    'OmegaMatter' : header['Omega0'],
                    'OmegaBaryon' : header['OmegaBaryon'],
                    'OmegaLambda' : header['OmegaLambda']}            

        numpartval   = header['NumPart_Total'][0]
        numpartunit  = self.inputfunc.set_unit(vardescription="Total number of gas particles in the simulation", 
                                Lunit=1, 
                                aFact=0, 
                                hFact=0)
        numpart      = {'Value':numpartval,'Info':numpartunit}
        numpart['Info'].pop('CGSConversionFactor',None)

        initMassTable = header['MassTable']

        DMmass        = initMassTable[1]
        OmegaDM       = (cosmo['OmegaMatter']-cosmo['OmegaBaryon'])
        Gasmass       = DMmass * (cosmo['OmegaBaryon']/OmegaDM)

        Dmunit        = self.inputfunc.set_unit(vardescription="Initial Dark Matter mass", 
                                Lunit=1.989e43,                                     
                                aFact=0, 
                                hFact=-1)
        Gmunit        = self.inputfunc.set_unit(vardescription="Initial Gass particle mass", 
                                Lunit=1.989e43,                                     
                                aFact=0, 
                                hFact=-1)
        DMmass        = {'Value':DMmass,'Info':Dmunit}
        DMmass  = self.inputfunc.assing_unit_unyt(DMmass,'Masses')
        DMmass       = {'Value':DMmass,'Info':Dmunit}
        DMmass['Info'].pop('CGSConversionFactor',None)


        Gasmass       = {'Value':Gasmass,'Info':Gmunit}
        Gasmass  = self.inputfunc.assing_unit_unyt(Gasmass,'Masses')
        Gasmass       = {'Value':Gasmass,'Info':Gmunit}
        Gasmass['Info'].pop('CGSConversionFactor',None)

        MassTable     = {'DarkMatterMass':DMmass,'GasMass':Gasmass}
                            
        # compute some extra variables
        H0            = h * 100 * 1e5 * unyt.cm / unyt.s / constants["Mpc"]       # H0 in 1/s
        rhoc          = 3*H0**2 / (8*np.pi*constants["Grav"]) # critical density in g/cm^3
        cosmo["H0"]   = H0
        cosmo["rhoc"] = rhoc
        cosmo["rhob"] = rhoc * cosmo["OmegaBaryon"]
            
        Header = {'BoxSize':box, 'Cosmo' : cosmo, 'NumPartTot':numpart, 'MassTable':MassTable}

        return Header

    def read_particles(self):
        
        groupname = self.groupname
        groupdic  = self.groupdic
        header    = self.read_header() 
        sightline = self.fileparams['sightline']
        
        sightline = self.inputfunc.set_sightlineparams(sightline,header)
        
        if self.snaptype == 'los':
            
            print("Not supported!")
        if self.snaptype == 'snapshot':
            
            axis_dic,xpos,ypos,zpos,los_length,msl_x3  = self.inputfunc.BuildLOS(sightline,header)
            sim_xmin,sim_xmax = axis_dic[0]
            sim_ymin,sim_ymax = axis_dic[1]
            sim_zmin,sim_zmax = axis_dic[2]
        
            # To read a region in hydrangea you have to define it as a rectangle, 
            # and you have to define where it is centered stp and it's  length given
            # by half of the sides.
            axis_2lenght = los_length * 0.5            
            stp = np.array([sim_xmin+msl_x3,sim_ymin+msl_x3,sim_zmin+axis_2lenght])

            if sightline['z-axis']==0:
                stp = np.array([sim_xmin+axis_2lenght,sim_ymin+msl_x3,sim_zmin+msl_x3])
                otp = np.array([axis_2lenght,msl_x3,msl_x3])             
            elif sightline['z-axis']==1:
                stp = np.array([sim_xmin+msl_x3,sim_ymin+axis_2lenght,sim_zmin+msl_x3])
                otp = np.array([msl_x3,axis_2lenght,msl_x3])
            else:
                otp = np.array([msl_x3,msl_x3,axis_2lenght])
                
            self.HY_snap = hy.ReadRegion(self.fname, part_type=0,anchor=stp,size= otp,shape='box',units='data')
            SmoothingL   = self.HY_snap.read_data('SmoothingLength') 
            Positions    = self.HY_snap.read_data("Coordinates")

            mask_x     =  ((Positions[:,sightline['x-axis']]>(xpos-SmoothingL)) & (Positions[:,sightline['x-axis']]<(xpos+SmoothingL)))
            mask_y     =  ((Positions[:,sightline['y-axis']]>(ypos-SmoothingL)) & (Positions[:,sightline['y-axis']]<((ypos+SmoothingL))))
            mask_z     =  ((Positions[:,sightline['z-axis']]>zpos) & (Positions[:,sightline['z-axis']]<((zpos+los_length))))

            self.los_mask = np.where((mask_x) & (mask_y) & (mask_z))[0]

            
            
        particles = {'Densities'        : self.read_variable(varname = groupname + '/' + groupdic['Densities']),
                     'SmoothingLengths' : self.read_variable(varname = groupname + '/' + groupdic['SmoothingLengths']),
                     'Masses'           : self.read_variable(varname = groupname + '/' + groupdic['Masses']),
                     'Positions'        : self.read_variable(varname = groupname + '/' + groupdic['Positions']),
                     'Velocities'       : self.read_variable(varname = groupname + '/' + groupdic['Velocities']),
                     'Temperatures'     : self.read_variable(varname = groupname + '/' + groupdic['Temperatures']),
                     'Elements'         : None, 
                     'Abundances'       : None  }
            
        elements,abundances,metallicity    = self.inputfunc.set_fractions_metallicity(self.read_variable,particles)    
        
        particles['Elements']   = elements
        particles['Abundances'] = abundances
        particles['Metallicities'] = metallicity
        
        if self.fileparams['ionparams']['SFR_properties']['modify_particle']:            
            particles['StarFormationRate']   = self.inputfunc.set_SFR(self.read_variable,groupname,particles)

        
        particles['Positions']['Value'] = self.inputfunc.flip_los(particles['Positions']['Value'],sightline)
        particles['Velocities']['Value'] = self.inputfunc.flip_los(particles['Velocities']['Value'],sightline)        
        
        return particles,sightline
        
        
    def read_variable(self,varname='PartType0/Density'):
        """
            Reads a particular hdf5 group and formats it's attributes

            Parameters
            ----------

            varname : str
            name of the hdf5 group with full hdf5 path
        """

        
        
        info   = self.inputfunc.read_group(groupname = varname)
        
        if self.snaptype == 'los':
            print("Not supported")
        
        elif self.snaptype == 'snapshot':
            varname_frmtd = varname.replace("PartType0/","")
            values = self.HY_snap.read_data(varname_frmtd)[self.los_mask]

        return {'Value': values, 'Info': info}

    
    
    

    
    
    
        
    
