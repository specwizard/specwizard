
import numpy as np
import h5py
from SpecWizard_BuildInput import Build_Input
from scipy.interpolate import interp1d
import SimulationInputKeys
import SimulationInputKeys as Skeys
import os
from SpecWizard_Elements import Elements
from SpecWizard_Input import ReadData
from SpecWizard_ProjectData import SightLineProjection
import Phys
constants = Phys.ReadPhys()
from SpecWizard_ComputeOpticaldepth import ComputeOpticaldepth
import random as rd
import copy
from SpecWizard_Lines import Lines

class LongSpectra:
    
    def __init__(self,wizard={}):
        
        self.wizard     = wizard
        self.lambda_min = wizard['longspectra']['lambda_min']
        self.lambda_max = wizard['longspectra']['lambda_max']
        self.dlambda    = wizard['longspectra']['dlambda']
        self.z_qsr      = wizard['longspectra']['z_qsr']
        self.delta_z    = wizard['longspectra']['delta_z']
        self.file_dir   = wizard['longspectra']['file_dir']
        self.pixkms     = 1
        self.npix       = int((self.lambda_max-self.lambda_min)/self.dlambda)
        self.wavelength =  self.lambda_min + (np.arange(self.npix) * self.dlambda) 

        
        self.paper = True

    def vel_from_wv(self,lambda_f,lambda_min,c_kms=1):
        '''
        Calculates velocity [km/s] from wavelength range
        '''
        return np.log(lambda_f/lambda_min) * c_kms

    def find_index(self,array2find,value):
        return int(np.argmin(abs(array2find-value)))

    def get_hdf5_files(self,los_dir):
        '''
        This function returns a list with all the hdf5 files in the directory assuming this are the files that will be use for building the longspectrum 
        '''
        los_files = []
        for file in os.listdir(los_dir):
            if file.endswith("hdf5"):
                los_files.append(file)
        return los_files
    
    def contaminant_lambda0(self):

        atomfile_name = self.wizard['ionparams']['atomfile']
        atomfile = h5py.File(atomfile_name,'r')
        zqsr = self.z_qsr
        l_min = self.lambda_min
        l_max = self.lambda_max
        lambda0z = {}
        # we load all the available line information that we have and we link their f-value,ion and element name. 
        for element in atomfile.keys():

            for ion in atomfile[element].keys():
                lambda0s = atomfile[element][ion]['lambda0'][...]
                fvalues  = atomfile[element][ion]['f-value'][...]
                for lambda0,fvalue in zip(lambda0s,fvalues):
                    lambda0z[lambda0] = {}
                    lambda0z[lambda0]['fvalue'] = fvalue
                    lambda0z[lambda0]['ion']    = ion
                    lambda0z[lambda0]['element'] = element


        all_lambda0 = np.array(list(lambda0z.keys()))
        # We want all the lines that at z_quasar lambda_0> lambda_min and lambda_0 < lambda_max and that at z=0 lambda_0<lambda_max

        self.lambda0_dic = lambda0z 
        mask1  = ((np.array(all_lambda0) * (1+zqsr)) > l_min)
        mask2  = ((np.array(all_lambda0) < l_max))
        total_mask = mask1 & mask2 

        contaminant_lambda0 = all_lambda0[total_mask]    
    
        self.contaminant_lambda0s = contaminant_lambda0
    
    def check_if_ion_contaminates(self,ions_we_want):
        ions2do = []
        cont_contrib_ions = np.unique([self.lambda0_dic[key]['ion'] for key in self.contaminant_lambda0s])

        for element,ion in ions_we_want:
            if ion in cont_contrib_ions:
                ions2do.append((element,ion))
            else:
                print(ion, "Does  not contribute in this wavelength range.")

        return ions2do
    def get_z_from_files(self,los_files):
        '''
        This function returns a dictionary containing the file name, and to what redshift this file corresponds
        '''
        aux_wizard = self.wizard.copy()
        sim_type = aux_wizard['file_type']['sim_type']
        sim_dics  = Skeys.get_simkeys((sim_type))
        los_dir   = aux_wizard['longspectra']['file_dir']
        file_los = {}
        for i,lfile in enumerate(los_files):

            file_los[i] = {}

            whereisZ = sim_dics['Header']['Redshift']
            lodfile = h5py.File(los_dir+lfile,'r')
            aux_wizard = copy.deepcopy(self.wizard)
            aux_wizard['file_type']['snap_type'] = 'los' 
            aux_wizard['snapshot_params']['directory'] = los_dir
            aux_wizard['snapshot_params']['file'] = lfile
            #
            file_los[i]['file'] = lfile
            file_los[i]['snap_type'] = 'los' 
            file_los[i]['los_dir']   = los_dir 
            file_los[i]['Redshift'] = lodfile[whereisZ[0]].attrs.get(whereisZ[1])    

        return file_los
    
    
    def Rebin(self,spectrum, wavelong,waveshort):
        ''' rebin the function y(x) to the new bins xnew
         Interpolation is performed such that the mean of the function is conserved
         Input: 
           -spectrum: dictionary, containing
              -L:       linear extent of the wavelength range
              -wave:    wavelength or Hubble velocity
              -flux:    flux
           -wave:       new wavelength or Hubble velocity to rebin to 
           '''
        # determine pixel size
        L    = waveshort.max() #spectrum["L"]
        npix = len(spectrum) #spectrum["npix"]
        pix  = L/npix
        x0   = np.arange(npix) * pix
        vHubble = waveshort
        dx = vHubble[1:] -vHubble[0:-1]
        dx = np.concatenate([dx, [L - x0[-1]]])
        #
        xR = np.copy(vHubble) + dx  # right side of pixel

        # cumulative sum of the function y0(x0)
        f = np.cumsum(spectrum * dx)

        # new pixels
        dxnew = wavelong[1:] - wavelong[0:-1]  # size of pixels
        dxnew = np.concatenate([dxnew, [L - wavelong[-1]]])
        xnewR = np.copy(wavelong) + dxnew  # right side of pixels

        # interpolation function
        finterp = interp1d(xR,
                           f,
                           kind='linear',
                           bounds_error=False,
                           fill_value=(0, f[-1]))
        fnew = finterp(xnewR)

        # rebinned value
        fbinned = np.concatenate([[0], fnew])
        ynew = (fbinned[1:] - fbinned[0:-1]) / dxnew
        #
        newspectrum         = spectrum.copy()
        #newspectrum["wave"] = np.copy(wave)
        #newspectrum["flux"] = ynew
        return ynew
    
    
    
    def get_file_and_los(self,z=0,random=True):
        '''
        Choose randomly a file and a sightline within the z threshold 
        '''
        # From the directory in the Wizard dictionary we get all the files. 
        hdf5files = self.get_hdf5_files(self.file_dir)
        # creates a dictionary that contains the file name and the correspondent redshift
        los_dict  = self.get_z_from_files(hdf5files)
        # get all the redshifts from the los files 

        los_z      = np.array([los_dict[i]['Redshift'] for  i in range(len(los_dict))])

        #We get all the files that are within our delta_z of our current z
        mask = np.where(np.abs(los_z-z) < self.delta_z)[0]

        # We pick a random file and a random sightline from it, check what to do if we are using snapshot

        if len(mask)==0:

            return None
        if random:
            random_file  = rd.choice(mask) 
            random_sight = rd.randint(0,99)

        else:
            random_file  = mask[0]
            random_sight = 1
        los_dict[random_file]['nsight'] = random_sight
        return los_dict[random_file]
    
    def adjust_wizard(self,los_dic):
        '''
        We adjust the original wizard dictionary to the specifics of the sightline and z
        '''
        BuildInput = Build_Input()

        snap_type = los_dic['snap_type']
        snapdir   = los_dic['los_dir']
        snapfile  = los_dic['file']
        nsight    = los_dic['nsight']
        BuildInput.FileType(sim_type='eagle',snap_type=snap_type)
        BuildInput.SnapshotParams(path=snapdir,file=snapfile)

        table_type     = self.wizard['ionparams']['table_type']
        iondir         = self.wizard['ionparams']['iondir']
        fname          = self.wizard['ionparams']['fname']
        SFR_properties = self.wizard['ionparams']['SFR_properties']
        ElementIons    = self.wizard['ionparams']['Ions']
        ElementNames   = np.array(ElementIons)[:,0]
        atomdat        = self.wizard['ionparams']['atomfile']
        transitions = Elements(atomdat)
        BuildInput.SetIonTableParams(table_type=table_type,iondir=iondir,ions=ElementIons,fname=fname,SFR_properties=SFR_properties,atomfile=atomdat)

        wizard    = BuildInput.Sightline(nsight=nsight)

        return wizard

    def insert_in_long_spectra(self,long_spectra,snapshot,projected_LOS,opticaldepth,ions2do,roll=True):
        
        pixel_kms = snapshot.ToCGS(projected_LOS["pixel_kms"]) / 1e5
        npix      = projected_LOS["npix"]
        vel_kms   = np.arange(npix) * pixel_kms 
        box_cgs   = snapshot.ToCGS(snapshot.header["BoxSize"])[0]       # in proper cm
        box_kms   = box_cgs * snapshot.Hubble() / 1e5 
        z_sim     = snapshot.header['Cosmo']['Redshift']

        lambda_min = self.lambda_min
        c_kms      = constants['c']/ 1e5
        amount_rolled  = 0
        for ion in ions2do:

            if self.paper== True:
                if self.readion == True:
                    if ion ==  ('Hydrogen', 'H I'):
                        opticaldepth[ion]['Optical depths'] = opticaldepth['SimIons'][ion]['Optical depths']

            tau         = snapshot.ToCGS(opticaldepth[ion]['Optical depths'])
            lambda0     = opticaldepth[ion]['lambda0']
            fvalue      = opticaldepth[ion]['f-value']
            #We use the periodic boundaries so we can extend 3 times the spectra, this to avoid problems with sharp edges
            extended_tau = np.concatenate((tau,tau))
            extended_tau = np.concatenate((tau,extended_tau))
            extended_tau  = extended_tau
            if roll:
                if amount_rolled ==0:
                    nindx = len(extended_tau)-1
                    amount_rolled = nindx-np.where(extended_tau==extended_tau.min())[0].max()
                    rolled = np.roll(extended_tau,amount_rolled)
                else:
                    rolled = np.roll(extended_tau,amount_rolled)
                    
                extended_tau = rolled
            
            
            llinelambda_start = lambda0 * (1 +z_sim) 
            llinelambda_end   =  llinelambda_start * np.exp(box_kms/c_kms)

            line_velstart =   self.vel_from_wv(llinelambda_start,lambda_min,c_kms)
            line_velend   =   self.vel_from_wv(llinelambda_end,lambda_min,c_kms)
            velarr_extended = np.linspace(line_velstart-npix,line_velend+npix,3*npix)

            rebinned   = self.Rebin(extended_tau,self.velocity_array,velarr_extended)
            indx_in    = self.find_index(self.velocity_array,line_velstart)
            indx_fn    = self.find_index(self.velocity_array,line_velend)
            
            if len(rebinned[indx_in:indx_fn]) == 0:
                #continue
                long_spectra['Ions'][ion]["OD"] += 0 # np.zeros_like(long_spectra['Ions'][ion]["OD"])
                
            else:
                long_spectra['Ions'][ion]["OD"][indx_in:indx_fn] += rebinned[indx_in:indx_fn]
            
            
            if long_spectra['Ions'][ion]["lambda0"] == 0:
                long_spectra['Ions'][ion]["lambda0"] = lambda0 
                long_spectra['Ions'][ion]["f-value"] = fvalue 
        return long_spectra

    def do_long_spectra(self):
        wizard = self.wizard
        delta_z = self.delta_z
        z_qsr   = self.z_qsr
        z       = 0
        c_kms      = constants['c'] / 1e5
        self.contaminant_lambda0()

        #We define the long spectra array in velocity space
        velocity_start = 0
        velocity_end   = self.vel_from_wv(self.lambda_max,self.lambda_min,c_kms)
        velocity_npix  = int((velocity_end-velocity_start) / self.pixkms)
        velocity_array = velocity_start + (np.arange(velocity_npix) * self.pixkms)
        #we create the optical depth array for the long spectra
        long_tau       = np.zeros_like(velocity_array,dtype=np.float)
        
        ions_we_want     = wizard['ionparams']['Ions']

        ions2do = self.check_if_ion_contaminates(ions_we_want)
        self.velocity_array = velocity_array
        long_spectra   = {}
        
        long_spectra['velocities'] = velocity_array
        long_spectra['Ions']       = {}


        for ions in ions2do:
            long_spectra['Ions'][ions] = {}
            long_spectra['Ions'][ions]["OD"] = long_tau.copy()
            long_spectra['Ions'][ions]["lambda0"] = 0
            long_spectra['Ions'][ions]["f-value"]= 0
        
         
        while z<z_qsr:
            print('Doing z: ',z)
            # We get choose randomly or not, a file and a los that is within dz 
            los_dict = self.get_file_and_los(z,random=True)

            #This means it did not found a file for the current z 
            if los_dict == None:
                print("No more files found within dz")
                z = z_qsr
                continue
        # we adjust the wizard dictionary to have the specifics of the choosen los
            wizard = self.adjust_wizard(los_dict)
            self.readion = False

            if self.paper==True:
                crit = abs(z - 3.017)
                print(crit)
                
                if crit <0.045:
                    self.paper = False
                    self.readion = True
                    print("this is happening")
                    wizard['file_type'] = {'sim_type': 'eagle', 'snap_type': 'snapshot'}
                    wizard['snapshot_params'] = {'directory': '/madfs/data/dc-syke1/Eagle/ScienceRuns/Planks1/L0050N1504/PE/RECALIBRATED/data/snapshot_012_z003p017/','file': 'snap_012_z003p017.0.hdf5'}
                    x_pos,y_pos,z_pos = 14.154414,9.485295 ,0 
                    wizard['sightline']['ProjectionAxes']= ['simx', 'simy', 'simz']
                    wizard['sightline']['x-axis']= 0
                    wizard['sightline']['y-axis']= 1
                    wizard['sightline']['z-axis']= 2
                    wizard['sightline']['ProjectionStart'] = (x_pos/50.0,y_pos/50.0,0)
                    wizard['sightline']['x-position'] = x_pos/50.0
                    wizard['sightline']['y-position'] = y_pos/50.0
                    wizard['sightline']['z-position'] = 0
                    wizard['sightline']['ProjectionLength']=1
                    wizard['sightline']['ProjectionExtend'] = {'extend': False, 'extendfactor': 3}
                    wizard['extra_parameters'][ 'ReadIonFrac'] = {'ReadIonFrac': True,
                                                                    'ReadHydrogen': True,
                                                                    'HI': 'HydrogenOneFraction',
                                                                    'ReadHelium': False,
                                                                    'He': '',
                                                                    'fname_urchin': '/madfs/data/dc-syke1/Urchin/Runs/ScienceRuns/Planck1/L0050N1504/RECAL/data/snapshot_012_z003p017/Combine/urchin_snap_012_z003p017.0.hdf5'}
                    print(wizard['sightline'])
            snapshot  = ReadData(wizard = wizard)
                
            particles = snapshot.read_particles()

            # We re-scale relevant parameters for the redshift we want, from here on z_sim = z

            snapshot.header['Cosmo']['Redshift'] = z 
            z_sim     = snapshot.header['Cosmo']['Redshift']
            box_sim   = snapshot.ToCGS(snapshot.header['BoxSize'])[0]
            #We calculate the extend in redshfit of the sightline
            dz_sim    = (1+z_sim) *  (np.exp(snapshot.Hubble() * box_sim  / constants['c'])-1)


            sightlineprojection  = SightLineProjection(wizard)
            projected_LOS        = sightlineprojection.ProjectData(particles)

            wizard['ODParams']['VoigtOff'] = True
            cspec                = ComputeOpticaldepth(wizard)
            opticaldepth         = cspec.MakeAllOpticaldepth(projected_LOS)

            long_spectra = self.insert_in_long_spectra(long_spectra,snapshot,projected_LOS,opticaldepth,ions2do,roll=True)        
            
            z += dz_sim

        return long_spectra
    
    def add_contaminants(self,longspectra):
        c_kms   = constants['c'] / 1e5
        wizard  = self.wizard
        atomdat = wizard['ionparams']['atomfile']
        transitions = Elements(atomdat)
        
        ions2do = (longspectra['Ions']).keys()
        velocity_array = longspectra['velocities']
        lambda0_dic = self.lambda0_dic
        contlambda = self.contaminant_lambda0s
        all_lines =  {}

        for (element,ion) in ions2do:
            if element not in all_lines.keys():
                all_lines[element] = {}

            linevals = [[lambda0,lambda0_dic[lambda0]['fvalue']]  for lambda0 in contlambda if lambda0_dic[lambda0]['ion']==ion]

            strongest_l0 = longspectra['Ions'][(element,ion)]["lambda0"]
            strongest_fv = longspectra['Ions'][(element,ion)]["f-value"]

            velocity_start  = self.vel_from_wv(strongest_l0,self.lambda_min,c_kms)
            velstart_indx   = self.find_index(velocity_array,velocity_start)

            last_indx      = np.where(longspectra['Ions'][ (element,ion)]["OD"]>0)[0][-1]
            len_arr        = last_indx-velstart_indx+1
            end_indx       = velstart_indx+len_arr

            tau_strong     = longspectra['Ions'][(element,ion)]["OD"][velstart_indx:end_indx] / (strongest_l0 * strongest_fv)
            
            for (lambda0,fvalue) in linevals:

                scale_tau       = tau_strong * (lambda0*fvalue) 
                temp_tau        = np.zeros_like(longspectra['Ions'][(element,ion)]["OD"])
                velocity_start  = self.vel_from_wv(lambda0,self.lambda_min,c_kms)
                scale_tau,velstart_indx,end_indx = self.velocity_index_check(scale_tau,velocity_start,velocity_array,len_arr)
                    
                temp_tau[velstart_indx:end_indx] = scale_tau

                all_lines[element][lambda0] = temp_tau
            
            if False:#(element,ion) ==('Hydrogen','H I'):
               all_lines = self.add_dw_and_lls(all_lines,velocity_array,len_arr)
            total_tau = np.zeros_like(temp_tau)
            
            for line in all_lines[element].keys():
                total_tau += all_lines[element][line]

            all_lines[element]['total'] = total_tau
        
        return_dic = {}
        return_dic["Ions"] = {}
        return_dic["Ions"] = all_lines
        return return_dic
        
    def velocity_index_check(self,scale_tau,velocity_start,velocity_array,len_arr):
        if velocity_start < 0:
            velstart_indx = 0 
            offindx       = int(abs(velocity_start))
            scale_tau     = scale_tau[offindx:]
            end_indx      = int(len(scale_tau))
            
        else:
            velstart_indx   = self.find_index(velocity_array,velocity_start)
            end_indx        = velstart_indx + len_arr
            
            if end_indx > len(velocity_array):
                end_indx = len(velocity_array)-1
                contr_indx = int(end_indx-velstart_indx)
                scale_tau = scale_tau[:contr_indx]
        return scale_tau,velstart_indx,end_indx

    def add_dw_and_lls(self,all_lines,velocity_array,len_arr):
        hi_phi = all_lines['Hydrogen'][1215.6701]
        temp_phi = np.zeros_like(hi_phi).copy()
        lines = Lines( v_kms =velocity_array ,box_kms=velocity_array.max())

        print("Calculating H I damping wings:...")
        dw      =  lines.convolvelorentz(hi_phi)
        print("Done!")

        ll_lambda0 = 912
        hi_lls  =  lines.convolveLymanLimit(hi_phi)
        print(len(hi_lls))
        c_kms = constants['c']/ 1e5
        velocity_start  = self.vel_from_wv(ll_lambda0,1215.6701,c_kms)
        print(velocity_start,self.lambda_min,c_kms)

        hi_lls,velstart_indx,end_indx = self.velocity_index_check(hi_lls,velocity_start,velocity_array,len_arr)
        print(len(hi_lls))
        print(velstart_indx,end_indx) 
        temp_phi[velstart_indx:end_indx] = hi_lls

        all_lines['Hydrogen'][1215.6701] = dw
        all_lines['Hydrogen']["Ly_limit_system"] = temp_phi

        return all_lines
        
    def rebin_to_spectrograph(self,long_spectra):
        c_kms = constants['c'] / 1e5
        velocity_array = self.velocity_array
        long_spectra = long_spectra["Ions"]
        wavelength_fine_array = self.lambda_min * np.exp(velocity_array/c_kms)
        wavelength_tau = {}
        lines  = list(long_spectra.keys())
        if "OD" in  long_spectra[lines[0]]:
            wavelength_tau['lines'] = {}

            for line in lines:
                wavelength_tau['lines'][line] = wave_rebin_total = self.Rebin(long_spectra[line]["OD"],self.wavelength,wavelength_fine_array)        
        else:
            for element in lines:
                wavelength_tau[element] = {}
                wavelength_tau[element]['lines'] = {}

                for line in long_spectra[element].keys():
                    wavelength_tau[element]['lines'][line] = wave_rebin_total = self.Rebin(long_spectra[element][line],self.wavelength,wavelength_fine_array)        

        wavelength_tau['wavelength'] = self.wavelength      

        return wavelength_tau