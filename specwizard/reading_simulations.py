import h5py 
import glob
import Phys
import re
import numpy as np
import swiftsimio as sw
import pyread_eagle as read_eagle
import hydrangea as hy
from SimulationInputKeys import get_simkeys
from SpecWizard_Elements import Elements

# physical constants in cgs units
constants = Phys.ReadPhys()

class InputFunctions:
    
    def __init__(self,fileparams={},header=None):
        '''
        This class containes a set of functions that are commonly shared among the different classes that read the data from simulations.
        
        input: Wizard, the dictionary created by SpecWizard_BuildInput
        '''
        
    
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
        
    def read_group(self, groupname ='Header'):
    # read all entries for this particular hdf5 group
        hfile = h5py.File(self.fname, "r")
        group    = hfile[groupname]
        grp_dict = {}
        for k in group.attrs.keys():
            grp_dict[k]= group.attrs[k]

        hfile.close()

        return dict(grp_dict)
    
    def set_unit(self, vardescription = 'text describing variable', Lunit=constants['Mpc'], aFact=1.0, hFact=1.0):
        '''
        "Formats the information for converting to CGS units" 
        '''
        return {'VarDescription': vardescription, 'CGSConversionFactor':Lunit, 'aexp-scale-exponent' :aFact, 'h-scale-exponent': hFact}

    
    def Hubble(self):
        ''' Hubble constant at redshift z, in units 1/s'''
        Hz  = self.header["Cosmo"]["H0"]
        z   = self.header["Cosmo"]["Redshift"]
        Hz *= np.sqrt(self.header["Cosmo"]["OmegaMatter"] * (1.+z)**3 + self.header["Cosmo"]["OmegaLambda"])
        return Hz
    
    def FormatTxt(self,text):
        '''
        Format text into lower case and eliminate spaces
        '''
        return text.lower().replace(" ","")
    
    def set_sightlineparams(self,sightline,header):
        '''
        Adds to the Wizard['sightline'] dictionary extra information about of the sightline, specially if a LOS file is used. 
        '''
        
        self.header = header
        if self.snaptype == 'los':
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
            xpos   = sightline_s[groupdic['x-position']] / box_size
            ypos   = sightline_s[groupdic['y-position']] / box_size
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

        # compute the velocity extent of the simulation volume in the sightline direction
        z_axis  = sightline["z-axis"]        
        box_cgs = self.ToCGS(header["BoxSize"])[z_axis]       # in proper cm
        box_kms = box_cgs * self.Hubble() / 1e5                    # in km/s
        unit    = self.set_unit(vardescription='Extent of simulation volume in direction of projection in terms of Hubble velocity',
                        Lunit  = 1e5,
                        aFact  = 0.,
                        hFact  = 0.)
        
        sightline["Boxkms"] = {'Value': box_kms, 'Info': unit}
        #
        box      = header["BoxSize"]["Value"][z_axis]
        unit     = self.set_unit(vardescription=header["BoxSize"]["Info"]["VarDescription"],
                        Lunit  = header["BoxSize"]["Info"]['CGSConversionFactor'],
                        aFact  = header["BoxSize"]["Info"]['aexp-scale-exponent'],
                        hFact  = header["BoxSize"]["Info"]['h-scale-exponent'])
        box      = {'Value':box, 'Info': unit}
        sightline["Box"] = box
        
        sightline['short-LOS'] = False
        if sightline['ProjectionLength'] < 1:
            sightline['short-LOS'] = True

        return sightline
        
    def BuildLOS(self,sightline,header):
        '''
        Build LOS from full snapshot, for eagle uses ReadEagle and for swift swiftsim io. We define a region and later an additional mask.
        ouput: Will calculate self.snap and self.los_mask. The former contains a region defined by the simulation reading functions e,g ReadEagle.
        and the latter is a more fine mask to only include particles that are within their smoothing length to sightline. 
        '''
        
        def MeanSeparationLenght(N,V):
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
        '''
         Reads and formats the ion fraction information from simulations that have support.
        '''
        
        groupdic         = self.groupdic
        userions         = self.fileparams['ionparams']['Ions']

        if self.simtype == 'colibre':

            # In the case of colibre the ion fraction is saved as a fraction of hydrogen content, while specwizard uses ion as a fraction of the element of the ion. 
            atomfile = self.fileparams['ionparams']['atomfile']
            elements = Elements(atomfile)
            elementnames = [name[0] for name in userions] 
            elements_info = elements.ElementParameters(ElementNames=elementnames)

            IonFracs = {}
            colfile    = h5py.File(self.fname,'r')
            col_ions   = colfile['SubgridScheme']['NamedColumns']['SpeciesFractions'][...].astype('str')
            Abundances = particles['Abundances']
            dens_cgs   = self.ToCGS(particles['Densities'])
            nH         = dens_cgs * Abundances['Hydrogen']['Value'] / constants['mH']
            
            for elements, ion in userions:
                mass_e   =  elements_info[elements]['Weight'] * constants['amu']
                ne       = dens_cgs * Abundances[elements]['Value'] / mass_e
                ion_formated = self.FormatTxt(ion)
                col_ions_formated = np.array([self.FormatTxt(col_ion) for col_ion in col_ions])
                ions_indx = np.where(col_ions_formated==ion_formated)[0]
                if len(ions_indx) ==0:
                    continue
                ionname   = str(col_ions[ions_indx][0])
                IonFrac =  read_variable(varname = groupname + '/' + groupdic['IonFractions'] + '/'+ ionname)
                IonFrac['Value'] = IonFrac['Value'] * (nH/ne)
                
                IonFracs[ion] = {'Value': IonFrac['Value'],'Info': IonFrac['Info']}
        return IonFracs    
    
    def set_fractions_metallicity(self,read_variable,particles):
        '''
        Some simulations do not contain metallicity or element fraction. This fuction attempts to read the simulation data, if it fails it set some default values
        Sets metallicities as zero 
        Sets elemen fraction of hydrogen and helium to primordial values. 
        '''
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
        '''
        This function will put the long projected axis as the third column in the output.
        '''
        
        temp = data.copy()
        temp[:,0] = data[:,sightinfo['x-axis']]
        temp[:,1] = data[:,sightinfo['y-axis']]
        temp[:,2] = data[:,sightinfo['z-axis']]

        return temp
    
    
    def ToCGS(self, variable):
        ''' 
        return simulations values for this variable in proper cgs units (no h)
        '''
        return variable["Value"] * self.CGSunit(variable)
    
    
    
    def CGSunit(self, variable):
        ''' 
        Use the information in the variable to compute the factor needed to convert
        simulation values to proper, h-free cgs units.
        This is of the form
        proper value = simulation value * CGSunit, where
        CGSunit = CGSConversionFactor * h**hexpo * a**aexpo
        '''
        #header     = self.read_group()
        #dependence on expansion factor
        ascale     = (1./(1+self.header["Cosmo"]["Redshift"]))**variable["Info"]["aexp-scale-exponent"]

        # dependence on hubble parameter
        hscale     = (self.header["Cosmo"]["HubbleParam"])**variable["Info"]["h-scale-exponent"]
        
        #
        return variable["Info"]["CGSConversionFactor"] * ascale * hscale
    
    
    def set_SFR(self,read_variable,groupname,particles):
        '''
        Attempts to read star formation rates, if not creates a zero array
        '''
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
            
        return particlesSFR
    
    
    
    
class ReadEagle:

    
    def __init__(self,fileparams={}):
        
        '''
        Class that contains functions to read particle data from the eagle simulation
        Input: Wizard dictionary, the output from SpecWizard_BuildInput
        '''
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
        '''
        Reads header from the eagle simulation.
        '''

        header = self.inputfunc.read_group()
 
        # the unit and h-dependence of Eagle is not stated; We assume it is in cMpc/h
        boxsize  = np.array([1., 1.0, 1.0 ]) * header['BoxSize']
        boxunit  = self.inputfunc.set_unit(vardescription="Extent of simulation volume", 
                                Lunit=constants['Mpc'], 
                                aFact=1.0, 
                                hFact=-1.0)
        box      = {'Value':boxsize, 'Info':boxunit}
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
        Gasmass       = {'Value':Gasmass,'Info':Gmunit}
        MassTable     = {'DarkMatterMass':DMmass,'GasMass':Gasmass}
        
        # compute some extra variables
        H0            = h * 100 * 1e5 / constants["Mpc"]      # H0 in 1/s
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
            groupname    = groupdic['groupname'].format(sightline['nsight'])

        if self.snaptype == 'snapshot':
            
            
            sightline
            axis_dic,xpos,ypos,zpos,los_length,msl_x3  = self.inputfunc.BuildLOS(sightline,header)
            sim_xmin,sim_xmax = axis_dic[0]
            sim_ymin,sim_ymax = axis_dic[1]
            sim_zmin,sim_zmax = axis_dic[2]
        
        
            self.RE_snap = read_eagle.EagleSnapshot(self.fname)
            self.RE_snap.select_region(sim_xmin,sim_xmax,sim_ymin,sim_ymax,sim_zmin,sim_zmax)            
            SmoothingL = self.RE_snap.read_dataset(0,'SmoothingLength')
            Positions  = self.RE_snap.read_dataset(0,'Coordinates')


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

        
        particles['Positions']['Value'] = self.inputfunc.flip_los(particles['Positions']['Value'],sightline)
        particles['Velocities']['Value'] = self.inputfunc.flip_los(particles['Velocities']['Value'],sightline)        
        
        return particles,sightline
        
        
    def read_variable(self,varname='PartType0/Density'):
        
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

            info   = self.inputfunc.read_group(groupname = varname)

        return {'Value': values, 'Info': info}
        
    def IsVelocity(self, name):
        # Check whether we are reading the velocity variable
        words    = name.split('/')
        Velocity = False
        for word in words:
                if (word == 'Velocity') or (word == 'Velocities'):
                    Velocity = True
        return Velocity


    
class ReadSwift:
    
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
        
        header    = self.inputfunc.read_group()
        cosmology = self.inputfunc.read_group('Cosmology')
        units     = self.inputfunc.read_group('Units')
        boxsize   = np.array(header['BoxSize'])
        boxunit   = self.inputfunc.set_unit(vardescription="Extent of simulation volume", Lunit=units['Unit length in cgs (U_L)'][0], aFact=1.0, hFact=0.0) # co-moving, no-h
        box       = {'Value':boxsize, 'Info':boxunit}
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
        Gasmass       = {'Value':Gasmass,'Info':Gmunit}
        MassTable     = {'DarkMatterMass':DMmass,'GasMass':Gasmass}
        
        
        # compute some extra variables
        H0            = h * 100 * 1e5 / constants["Mpc"]      # H0 in 1/s
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
            load_region = [[sim_xmin*oneMpc,sim_xmax*oneMpc],[sim_ymin*oneMpc,sim_ymax*oneMpc],[sim_zmin*oneMpc,sim_zmax*oneMpc]]
            swif_mask.constrain_spatial(load_region)

            self.SW_snap = sw.load(self.fname, mask=swif_mask)
            SmoothingL = self.SW_snap.gas.smoothing_lengths.value
            Positions  = self.SW_snap.gas.coordinates.value


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
        
        
        FWHM = 0.362
        particles['SmoothingLengths']["Value"] /= FWHM
        print("We divide Swift's smoothing length by {0:1.3f} to convert from FWHM to extent of finite support".format(FWHM))

        if (self.simtype == 'colibre'):
            
            particles['SimulationIonFractions'] = self.inputfunc.ReadAndShapeIonFrac(self.read_variable,particles,groupname)
        
        if self.fileparams['ionparams']['SFR_properties']['modify_particle']:
            
            particles['StarFormationRate']   = self.inputfunc.set_SFR(self.read_variable,groupname,particles)

        particles['Positions']['Value'] = self.inputfunc.flip_los(particles['Positions']['Value'],sightline)
        particles['Velocities']['Value'] = self.inputfunc.flip_los(particles['Velocities']['Value'],sightline)        
        
        return particles,sightline
    
    
    def read_variable(self,varname='PartType0/Density'):
        
        if self.snaptype == 'los':
            hfile  = h5py.File(self.fname, "r")
            values = np.array(hfile[varname][...])
            hfile.close()
        
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
        varname_formated = re.sub(r"([A-Z])", r"_\1", varname)

        return varname_formated[1:].lower()

    
class ReadHydrangea:
    
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
        Gasmass       = {'Value':Gasmass,'Info':Gmunit}
        MassTable     = {'DarkMatterMass':DMmass,'GasMass':Gasmass}
                            
        # compute some extra variables
        H0            = h * 100 * 1e5 / constants["Mpc"]      # H0 in 1/s
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
        info   = self.inputfunc.read_group(groupname = varname)
        
        if self.snaptype == 'los':
            print("Not supported")
        
        elif self.snaptype == 'snapshot':
            varname_frmtd = varname.replace("PartType0/","")
            values = self.HY_snap.read_data(varname_frmtd)[self.los_mask]

        return {'Value': values, 'Info': info}

    
    
    

    
    
    
        
    
