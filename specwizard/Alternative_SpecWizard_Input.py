
import numpy as np
import Phys
import importlib
import pyread_eagle as read_eagle
import swiftsimio as sw
import hydrangea as hy
from SpecWizard_Elements import Elements
import h5py
Phys = importlib.reload(Phys)
from SimulationInputKeys import get_simkeys

# physical constants in cgs units
constants  = Phys.ReadPhys()



# Eagle runs
Eagle = True
if Eagle:
    filetype = {'simtype':'Eagle', 'snaptype':'LOS'}  

    # location of snapshot
    snapdir = '/cosma7/data/Eagle/ScienceRuns/Planck1/L0050N1504/PE/RECALIBRATED/data/los/'

    # name of snapshot file
    snapfile = 'part_los_z3.275.hdf5'

else:
    filetype = {'simtype':'Swift', 'snaptype':'LOS'}  

    # location of snapshot
    snapdir = '/cosma6/data/dp004/dc-elbe1/qla_nu/runs/QLANU_M060_L020N0256/'

    # name of snapshot file
    snapfile = 'los_0019.hdf5'


fileparams = {'file_type'   : {'sim_type' : 'Eagle', 'snap_type' : 'LOS'},
              'snapshot_params':{'directory':snapdir, 'file':snapfile}}


# sightline parameters: 
sightline  = {'x-axis': 0, 'y-axis':1, 'z-axis': 2, 'x-position':0.0, 'y-position':0.0, 'nsight' : 0}

class ReadSnap:
    '''
    Read SPH particles from Eagle (Gadget), Hydrangea(C-Eagle) or Swift output
    input: Wizardparamas(output from running SpecWizard_BuildInput)
    '''
    def __init__(self, fileparams=fileparams):
        #
        self.fileparams = fileparams
        self.fname      = fileparams["snapshot_params"]["directory"] + '/' + fileparams["snapshot_params"]["file"]
        self.simtype    = fileparams["file_type"]["sim_type"]
        self.snaptype   = fileparams["file_type"]["snap_type"]
        self.readIonFrac = fileparams['extra_parameters']['ReadIonFrac']['ReadIonFrac'] 
        
        # read default simulation parameters 
        self.header   = self.ReadHeader()
    
        self.fileparams["Header"] = self.header
    def ReadHeader(self):
        ''' Read header and cosmology information from file '''
        
        # store default data
        if self.simtype == 'eagle':
            # read header information and store in default format
            header = self.ReadGroup()
        
            # the unit and h-dependence of Eagle is not stated; We assume it is in cMpc/h
            boxsize  = np.array([1., 1.0, 1.0 ]) * header['BoxSize']
            boxunit  = self.SetUnit(vardescription="Extent of simulation volume", 
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
            numpartunit  = self.SetUnit(vardescription="Total number of gas particles in the simulation", 
                                    Lunit=1, 
                                    aFact=0, 
                                    hFact=0)
            numpart      = {'Value':numpartval,'Info':numpartunit}
            
        if (self.simtype == 'swift') or (self.simtype == 'colibre'):
            header    = self.ReadGroup()
            cosmology = self.ReadGroup('Cosmology')
            units     = self.ReadGroup('Units')
            boxsize   = np.array(header['BoxSize'])
            boxunit   = self.SetUnit(vardescription="Extent of simulation volume", Lunit=units['Unit length in cgs (U_L)'][0], aFact=1.0, hFact=0.0) # co-moving, no-h
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
            numpartunit  = self.SetUnit(vardescription="Total number of gas particles in the simulation", 
                                    Lunit=1, 
                                    aFact=0, 
                                    hFact=0)
            numpart      = {'Value':numpartval,'Info':numpartunit}          
        if self.simtype == 'hydrangea':
            # read header information and store in default format
            header = self.ReadGroup()
        
            # the unit and h-dependence of Eagle is not stated; We assume it is in cMpc/h
            boxsize  = np.array([1., 1.0, 1.0 ]) * header['BoxSize']
            boxunit  = self.SetUnit(vardescription="Extent of simulation volume", 
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
            numpartunit  = self.SetUnit(vardescription="Total number of gas particles in the simulation", 
                                    Lunit=1, 
                                    aFact=0, 
                                    hFact=0)
            numpart      = {'Value':numpartval,'Info':numpartunit}          
                            
        # compute some extra variables
        H0            = h * 100 * 1e5 / constants["Mpc"]      # H0 in 1/s
        rhoc          = 3*H0**2 / (8*np.pi*constants["Grav"]) # critical density in g/cm^3
        cosmo["H0"]   = H0
        cosmo["rhoc"] = rhoc
        cosmo["rhob"] = rhoc * cosmo["OmegaBaryon"]
            
        Header = {'BoxSize':box, 'Cosmo' : cosmo, 'NumPartTot':numpart}
        
        return Header
    
    def ReadParticles(self, specparams = None):
        
        
        ''' Read particles along a given sight line
        
        input: Wizardparams (output from Specwizard_BuildInput)
        sightline  = {'x-axis': 0, 'y-axis':1, 'z-axis': 2, 'x-position':0.0, 'y-position':0.0, 'nsight' : 0}
        
        For LOS output we will simply read sight line number nsight from the file
        
        If snaptype is not LOS, we read a full snapshot. Which sight line to extract is set as follows:
        
        The sightline is parallel to the z-axis, and goes through the pt with (x,y) coordinates (x-position, y-position)
        (with units the same as those of the simulation)
        The integers x-axis, y-axis, z-axis decide which of the coordinates of the particle, correspond to this new set of axes.

        For example, for (x-axis=2, y-axis=1, z-axis=0), the sight line will - in the coordinates of the simulation,
        be along the x-axis, and go through the point with (y,z) = (x-position, y-position)
        
        
        Output:
           This routine outputs a dictionary
           'Particles' are the proerties of the SPH particles read in
           'SightInfo' is extra information on this sight line
        
        '''
        sightline = specparams['sightline']
        sightline['short-LOS'] = False
        sim_keys  = get_simkeys(self.simtype)
        groupdic  = sim_keys[self.snaptype]
        self.groupdic = groupdic
        # If we are using a los output we read the spacial information regarding the sightline.
        if self.snaptype == 'los':
            groupname    = groupdic['groupname'].format(sightline['nsight'])
            sightline_s  = self.ReadGroup(groupname = groupname)
            try:
                num_particles = sightline_s[groupdic['Number_of_part_this_los']]
            except: 
                num_particles = 0
            box_size = self.header["BoxSize"]['Value'][0]
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
            if self.simtype == 'swift':
                sightline['z-position'] = 0.
                sightline['x-axis']     = xaxis[0]
                sightline['y-axis']     = yaxis[0]
                sightline['z-axis']     = zaxis[0]  
            sightline['short-LOS'] = False
            sightline['ProjectionLength'] = 1

        #if snapshot output is used then we will construct the LOS using the BuildLOS function. 
        else:
            groupname  = groupdic['groupname']
            self.BuildLOS(sightline)



        particles = {'Densities'        : self.ReadVariable(varname = groupname + '/' + groupdic['Densities']),
                     'SmoothingLengths' : self.ReadVariable(varname = groupname + '/' + groupdic['SmoothingLengths']),
                     'Masses'           : self.ReadVariable(varname = groupname + '/' + groupdic['Masses']),
                     'Positions'        : self.ReadVariable(varname = groupname + '/' + groupdic['Positions']),
                     'Velocities'       : self.ReadVariable(varname = groupname + '/' + groupdic['Velocities']),
                     'Temperatures'     : self.ReadVariable(varname = groupname + '/' + groupdic['Temperatures']),
                     'Elements'         : None, 
                     'Abundances'       : None  }

        #We read the element fractions to do, if not found primordial quantities are used. 
        abundances = {}
        unit       = self.SetUnit(vardescription = 'mass fraction of element', 
                                            Lunit=1.0, aFact=0.0, hFact=0.0)

        elementnames     = groupdic['elementnames']
        userions         = self.fileparams['ionparams']['Ions']
        elements_we_want = np.array([userions[i][0] for i in range(len(userions))])
        elements2do      = np.intersect1d(elements_we_want,elementnames)
        
        try:    
            if len(elements2do) == 0:
                print("ERROR! No valid element or element name!")

            for elementname in elements2do:

                values  = self.ReadVariable(varname = groupname +'/'+ groupdic['ElementAbundance']+'/'+ elementname)['Value']
                info    = unit
                abundances[elementname] = {'Value':values, 'Info': unit}
            particles['Elements']   = elements2do
            particles['Abundances'] = abundances
        except:
            print("Element fraction not found using primordial quantities...")
            hydrogen    = np.zeros_like(particles['Densities']['Value']) + 1. - constants['Y']
            abundances['Hydrogen'] = {'Value':hydrogen, 'Info':unit}
            helium      = np.zeros_like(particles['Densities']['Value']) + constants['Y']
            abundances['Helium']  = {'Value':helium, 'Info':unit}
            particles['Elements']   = elements2do
            particles['Abundances'] = abundances

        #Check if metallicities are available. 
        try: 
            metallicities  = self.ReadVariable(varname = groupname + '/' + groupdic['Metallicities'])
            particles['Metallicities'] = metallicities
        except:
            print("Warning! Metallicities not found. Setting them to zero.")
            metallicity = np.zeros_like(particles['Densities']['Value'])
            unit        = self.SetUnit(vardescription = 'mass fraction of metals', 
                                       Lunit=1.0, aFact=0, hFact=0)
            particles['Metallicities'] = {'Value':metallicity, 'Info':unit} 

        #Will modify the star formation rate particles if requested

        if specparams['ionparams']['SFR_properties']['modify_particle']:
            try:
                
                SFR    = self.ReadVariable(varname = groupname + '/' + groupdic['StarFormationRate'])
                particles['StarFormationRate'] = SFR
            except:
                print('Warning! Not able to read Star formation Rate properties...setting SFR to zero!')
                SFR    = np.zeros_like(particles['Masses']['Value']) 
                unit        = self.SetUnit(vardescription = 'Star Formation Rate', 
                               Lunit=1.0, aFact=0, hFact=0)
                particles['StarFormationRate'] = {'Value':SFR, 'Info':unit}   

        if (self.simtype == 'swift') or (self.simtype == 'colibre'):        
            FWHM = 0.362
            particles['SmoothingLengths']["Value"] /= FWHM
            print("We divide Swift's smoothing length by {0:1.3f} to convert from FWHM to extent of finite support".format(FWHM))
#super uglY need to imporve!!


        if (self.simtype == 'colibre'):
            
            particles['SimulationIonFractions'] = self.ReadAndShapeIonFrac(particles,userions,groupname,groupdic)
        
        
        particles['Positions']['Value'] = self.flip_los(particles['Positions']['Value'],sightline)
        particles['Velocities']['Value'] = self.flip_los(particles['Velocities']['Value'],sightline)


        # compute the velocity extent of the simulation volume in the sightline direction
        z_axis  = sightline["z-axis"]        
        box_cgs = self.ToCGS(self.header["BoxSize"])[z_axis]       # in proper cm
        box_kms = box_cgs * self.Hubble() / 1e5                    # in km/s
        unit    = self.SetUnit(vardescription='Extent of simulation volume in direction of projection in terms of Hubble velocity',
                        Lunit  = 1e5,
                        aFact  = 0.,
                        hFact  = 0.)
        
        sightline["Boxkms"] = {'Value': box_kms, 'Info': unit}
        #
        box      = self.header["BoxSize"]["Value"][z_axis]
        unit     = self.SetUnit(vardescription=self.header["BoxSize"]["Info"]["VarDescription"],
                        Lunit  = self.header["BoxSize"]["Info"]['CGSConversionFactor'],
                        aFact  = self.header["BoxSize"]["Info"]['aexp-scale-exponent'],
                        hFact  = self.header["BoxSize"]["Info"]['h-scale-exponent'])
        box      = {'Value':box, 'Info': unit}
        sightline["Box"] = box
        
        self.fileparams['sightline'] = sightline
        return {'SightInfo': sightline, 'Particles' : particles, 'Header': self.header}
                
    def Hubble(self):
        ''' Hubble constant at redshift z, in units 1/s'''
        Hz  = self.header["Cosmo"]["H0"]
        z   = self.header["Cosmo"]["Redshift"]
        Hz *= np.sqrt(self.header["Cosmo"]["OmegaMatter"] * (1.+z)**3 + self.header["Cosmo"]["OmegaLambda"])
        return Hz
    
    def FormatTxt(self,text):
        return text.lower().replace(" ","")
        
    def SetUnit(self, vardescription = 'text describing variable', Lunit=constants['Mpc'], aFact=1.0, hFact=1.0):
        return {'VarDescription': vardescription, 'CGSConversionFactor':Lunit, 'aexp-scale-exponent' :aFact, 'h-scale-exponent': hFact}
  
    def ReadAndShapeIonFrac(self,particles,userions,groupname,groupdic):
        
        if self.simtype == 'colibre':
            
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
                IonFrac = self.ReadVariable(varname = groupname + '/' + groupdic['IonFractions'] + '/'+ ionname)
                IonFrac['Value'] = IonFrac['Value'] * (nH/ne)
                
                IonFracs[ion] = {'Value': IonFrac['Value'],'Info': IonFrac['Info']}
            return IonFracs
        
    def ReadVariable(self, varname = 'LOS1/Density'):
        if self.simtype == 'eagle':
            info   = self.ReadGroup(groupname = varname)
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
                
            elif self.snaptype == 'snapshot':
                varname_frmtd = varname.replace("PartType0/","")
                values = self.RE_snap.read_dataset(0,varname_frmtd)[self.los_mask]
                
        if (self.simtype == 'swift') or (self.simtype =='colibre'):
            check_for_elmt = varname.split('/')
            
            if check_for_elmt[1]==self.groupdic['ElementAbundance']:
                varname_aux = 'PartType0/ElementMassFractions'
                info_s  = self.ReadGroup(groupname = varname_aux)
                
            elif check_for_elmt[1]==self.groupdic['Metallicities']:
                varname_aux = 'PartType0/MetalMassFractions'
                info_s  = self.ReadGroup(groupname = varname_aux)

            elif check_for_elmt[1]==self.groupdic['StarFormationRate']:
                varname_aux = 'PartType0/StarFormationRates'
                info_s  = self.ReadGroup(groupname = varname_aux)
            
            
            elif self.simtype == 'colibre':
                if check_for_elmt[1]==self.groupdic['IonFractions']:
                    varname_aux = 'PartType0/SpeciesFractions'
                    info_s  = self.ReadGroup(groupname = varname_aux)
                else:
                    info_s  = self.ReadGroup(groupname = varname)

            else:
                info_s  = self.ReadGroup(groupname = varname)
                    
                    
            info    = {'VarDescription': info_s['Description'] + info_s['Expression for physical CGS units'],
                       'CGSConversionFactor': info_s['Conversion factor to CGS (not including cosmological corrections)'][0],
                       'h-scale-exponent': info_s['h-scale exponent'][0],
                       'aexp-scale-exponent': info_s['a-scale exponent'][0]}
                
            if self.snaptype == 'los':
                hfile  = h5py.File(self.fname, "r")
                values = np.array(hfile[varname][...])
                hfile.close()
            
            elif self.snaptype == 'snapshot':
                varname_frmtd = varname.replace("PartType0/","")
                varname_frmtd = varname_frmtd.split("/")
                if varname_frmtd[0] == 'SmoothingLengths':
                    varname_frmtd = 'smoothing_lengths'
                    values = (getattr(self.SW_snap.gas,varname_frmtd).value)[self.los_mask]

                elif varname_frmtd[0] == self.groupdic['ElementAbundance']:
                    varname_frmtd[0] = 'element_mass_fractions'
                    varname_frmtd[1] =varname_frmtd[1].lower()
                    values = (getattr(getattr(self.SW_snap.gas,varname_frmtd[0]),varname_frmtd[1]).value)[self.los_mask]
                
                
                elif self.simtype == 'colibre':
                    if varname_frmtd[0] == self.groupdic['IonFractions']:
                        values = (getattr(getattr(self.SW_snap.gas,'species_fractions'),varname_frmtd[1]).value)[self.los_mask]

                    else:
                        varname_frmtd = varname_frmtd[0].lower()
                        values = (getattr(self.SW_snap.gas,varname_frmtd).value)[self.los_mask]

                else:
                    varname_frmtd = varname_frmtd[0].lower()
                    values = (getattr(self.SW_snap.gas,varname_frmtd).value)[self.los_mask]
                    
        if self.simtype == 'hydrangea':
            info   = self.ReadGroup(groupname = varname)
            if self.snaptype == 'los':
                hfile  = h5py.File(self.fname, "r")
                values = np.array(hfile[varname][...])
                hfile.close()
            elif self.snaptype == 'snapshot':
                varname_frmtd = varname.replace("PartType0/","")
                values = self.HY_snap.read_data(varname_frmtd)[self.los_mask]

                
        return {'Value': values, 'Info': info}
        
    def CGSunit(self, variable):
        ''' 
        Use the information in the variable to compute the factor needed to convert
        simulation values to proper, h-free cgs units.
        This is of the form
        proper value = simulation value * CGSunit, where
        CGSunit = CGSConversionFactor * h**hexpo * a**aexpo
        '''
        #dependence on expansion factor
        ascale     = (1./(1+self.header["Cosmo"]["Redshift"]))**variable["Info"]["aexp-scale-exponent"]

        # dependence on hubble parameter
        hscale     = (self.header["Cosmo"]["HubbleParam"])**variable["Info"]["h-scale-exponent"]
        
        #
        return variable["Info"]["CGSConversionFactor"] * ascale * hscale
        
        
    def ToCGS(self, variable):
        ''' 
        return simulations values for this variable in proper cgs units (no h)
        '''
        return variable["Value"] * self.CGSunit(variable)
    
    
    def flip_los(self,data,sightinfo):
        '''
        This function will put the long projected axis as the third column in the output.
        '''
        
        temp = data.copy()
        temp[:,0] = data[:,sightinfo['x-axis']]
        temp[:,1] = data[:,sightinfo['y-axis']]
        temp[:,2] = data[:,sightinfo['z-axis']]

        return temp 
    def BuildLOS(self,sightline):
        '''
        Build LOS from full snapshot, for eagle uses ReadEagle and for swift swiftsim io. We define a region and later an additional mask.
        ouput: Will calculate self.snap and self.los_mask. The former contains a region defined by the simulation reading functions e,g ReadEagle.
        and the latter is a more fine mask to only include particles that are within their smoothing length to sightline. 
        '''
        
        def MeanSeparationLenght(N,V):
            ninv = V/N
            return pow(ninv,1/3)
                    
        # We calculate the mean separation length to define the region that will be masked by the simulation reading function. 
        N = self.header['NumPartTot']['Value']
        V = np.prod(self.header['BoxSize']['Value'])
        msl_x3 = 3*MeanSeparationLenght(N,V)
        # We get the los_lenght as a fraction and multiply by the box size so we can extract particles from 
        # z-position to z-position+los_length 
        
        if sightline['ProjectionLength'] < 1:
            sightline['short-LOS'] = True

        los_length  = sightline['ProjectionLength']        
        xpos        = sightline['x-position']
        ypos        = sightline['y-position']
        zpos        = sightline['z-position']
        
        box         = self.header['BoxSize']['Value'][0]
    
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
        

        # The functions do not return a very narrow region, so we use the smoothing lengths to mask further the LOS.
        if self.simtype == 'eagle':
            
            self.RE_snap = read_eagle.EagleSnapshot(self.fname)
            self.RE_snap.select_region(sim_xmin,sim_xmax,sim_ymin,sim_ymax,sim_zmin,sim_zmax)            
            SmoothingL = self.RE_snap.read_dataset(0,'SmoothingLength')
            Positions  = self.RE_snap.read_dataset(0,'Coordinates')

    
        elif (self.simtype == 'swift') or (self.simtype == 'colibre'):

            swif_mask = sw.mask(self.fname)
            Sboxsize = swif_mask.metadata.boxsize[0]
            oneMpc =  Sboxsize/Sboxsize.value                   #Since swiftsimIO uses units we have to 
            load_region = [[sim_xmin*oneMpc,sim_xmax*oneMpc],[sim_ymin*oneMpc,sim_ymax*oneMpc],[sim_zmin*oneMpc,sim_zmax*oneMpc]]
            swif_mask.constrain_spatial(load_region)

            self.SW_snap = sw.load(self.fname, mask=swif_mask)
            SmoothingL = self.SW_snap.gas.smoothing_lengths.value
            Positions  = self.SW_snap.gas.coordinates.value
                        
        elif self.simtype == 'hydrangea':
            # To read a region in hydrangea you have to define it as a rectangle, 
            # and you have to define where it is centered stp and it's  length given
            # by half of the sides.
            axis_2lenght = los_length * 0.5            
            stp = np.array([sim_xmin+msl_x3,sim_ymin+msl_x3,sim_zmin+axis_2lenght])

            if sightline['z-axis']==0:
                otp = np.array([axis_2lenght,msl_x3,msl_x3])             
            elif sightline['z-axis']==1:
                otp = np.array([msl_x3,axis_2lenght,msl_x3])
            else:
                otp = np.array([msl_x3,msl_x3,axis_2lenght])
                
            self.HY_snap = hy.ReadRegion(self.fname, part_type=0,anchor=stp,size= otp,shape='box',units='data')
            SmoothingL   =  self.HY_snap.read_data('SmoothingLength') 

            Positions    = self.HY_snap.read_data("Coordinates")
        
        self.los_mask = np.where((Positions[:,sightline['x-axis']]>(xpos-SmoothingL)) & (Positions[:,sightline['x-axis']]<(xpos+SmoothingL))
                               & (Positions[:,sightline['y-axis']]>(ypos-SmoothingL)) & (Positions[:,sightline['y-axis']]<((ypos+SmoothingL)))
                               & (Positions[:,sightline['z-axis']]>zpos) & (Positions[:,sightline['z-axis']]<((zpos+los_length))))[0]


            
    def ReadGroup(self, groupname ='Header'):
        # read all entries for this particular hdf5 group
        hfile = h5py.File(self.fname, "r")
        group    = hfile[groupname]
        grp_dict = {}
        for k in group.attrs.keys():
            grp_dict[k]= group.attrs[k]
        
        hfile.close()
    
        return dict(grp_dict)
    
    def IsVelocity(self, name):
        # Check whether we are reading the velocity variable
        words    = name.split('/')
        Velocity = False
        for word in words:
                if (word == 'Velocity') or (word == 'Velocities'):
                    Velocity = True
        return Velocity    