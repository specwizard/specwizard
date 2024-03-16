
import numpy as np
from os.path import exists
from os import remove
import copy
import h5py
#
#from SpecWizard_Elements import Elements
#from SpecWizard_Input import ReadSnap
#from SpecWizard_SplineInterpolation import ColumnTable
#from SpecWizard_IonizationBalance import IonizationBalance
#from SpecWizard_ProjectData import SightLineProjection

    
class OpticalDepth_IO:
    ''' Methods to write and read files containing spectra and optical dept '''
    def __init__(self, wizard=None, create=True):
        self.dirname    = wizard["Output"]["directory"]
        self.fname      = wizard["Output"]["fname"]
        self.wizard = wizard
        self.nspec      = 0
        self.header     = wizard["Header"]
        
        #
        self.ions       = wizard["ionparams"]["Ions"]
        element_names   = np.array([ self.ions[i][0] for i in range(len(self.ions))])
        self.elements   = element_names
        
        # if create: create file and add header
        if create:
            file = self.dirname + self.fname
            if exists(file):
                print("Removing file ", file)
                remove(file)
            # write header
            self.WriteHeader()
        
            
    def ReadGroup(self, groupname='Header'):
        # read all entries for this particular hdf5 group
        hfile = h5py.File(self.dirname + self.fname, "r")
        group = hfile[groupname]
        grp_dict = {}
        for k in group.attrs.keys():
            grp_dict[k] = group.attrs[k]

        hfile.close()
        return dict(grp_dict)

    def WriteGroup(self, dictionary, groupname='Header'):
        ''' write all entries of this dictionary to groupname '''
        hfile = h5py.File(self.dirname + self.fname, "a")
        group = hfile.create_group(groupname)
        #
        for key, value in dictionary:
            group.attrs[key] = value

        hfile.close()

    def WriteVariable(self, variable = 1, varname='LOS1/Density'):
        # save all the attributes of a variable and its value
        # the input variable is a dicttionary of format
        # {"Value":values, "Info": info}
        # where - values is an np array
        #       - info is the dictionary 
        #          info = {'Vardescription':'string', 'CGSConvertsionFactor':1, 'h-scale-exponent':1, 'aexp-scale-exponent':1}

        # open file
        hfile  = h5py.File(self.dirname + self.fname, "a")
        
        # create and write dataset
        Values = variable["Value"]
        hfile.create_dataset(varname, data=Values)
        dset   = hfile[varname]

        # add attributes
        Info   = variable['Info']
        for (key, value) in Info.items():
            dset.attrs[key] = value
    
        hfile.close()

    def WriteHeader(self):
        ''' Write header and cosmology information '''
        # create file and store Header data
        file   = self.dirname + self.fname
        header = self.header
        #
        hfile = h5py.File(file, "a")
        group = hfile.create_group('Header')
        hfile.close()
        self.WriteVariable(varname = 'Header/Box', variable = header["BoxSize"])
        
        #
        hfile = h5py.File(self.dirname + self.fname, "a")
        for key, value in header["Cosmo"].items():
            hfile["Header"].attrs[key] = value
            
        # add link to original snapshot
        hfile["Header"].attrs["snapshot_directory"]   = self.wizard["snapshot_params"]["directory"]
        hfile["Header"].attrs["snapshot_file"]        = self.wizard["snapshot_params"]["file"]
        
        # table for computing ionization balance
        hfile["Header"].attrs["ionization_directory"] = self.wizard["ionparams"]["iondir"]
        
        # number of sightlines
        hfile["Header"].attrs["number_of_sightlines"] = 0

        hfile.close()
        
    def write_to_file(self, projections):
        ''' Add contents of this optical depth data to file     '''

        
        #
        nsight       = projections["nsight"]
        info         = self.wizard['sightline']
        projection   = projections["Projection"]
        opticaldepth = projections["OpticaldepthWeighted"]
        
        # update the number of sightlines done
        self.nspec += 1
        
        # open hdf5 file
        hfile = h5py.File(self.dirname + self.fname, "a")
        hfile['Header'].attrs.modify('number_of_sightlines',self.nspec)


        # create group for this line of sight
        groupname = 'LOS_{}'.format(nsight)
        try:
            hfile.require_group(groupname)
        except:
            print("error: group already exists")
            hfile.close()
            return

        # add attributes
        for key, value in info.items():
            try:
                hfile[groupname].attrs[key] = value
            except:
                continue
        Box = groupname + '/Box_kms'
        self.WriteVariable(info["Boxkms"], varname = Box)
        
        #
        variables    = ['Velocities', 'Densities', 'Temperatures']
        # create group for each element
        for element in self.elements:
            elementgroup = groupname + '/' + element
            hfile.require_group(elementgroup)
            
            # create group for element weighted properties
            group = elementgroup + '/Element-weighted'
            
            if not group in hfile:
                hfile.require_group(group)

                # add attributes
                for(key, value) in projection.items():
                    try:
                        hfile[group].attrs[key] = value
                    except:
                        continue
                    pix_kms = group + '/' + 'pixel_kms'
                    self.WriteVariable(projection["pixel_kms"], varname = pix_kms)

                # add element-weighted variables
                for variable in variables:
                    self.WriteVariable(projection["Element-weighted"][element][variable], 
                                       group+'/' + variable)
                
        # create ion group for each ion
        for transition in self.ions:
            (element, ion) = transition
            group = groupname + '/' + element + '/' + ion
            hfile.require_group(group)
            
           # ion-weighted properties
            group = groupname + '/' + element + '/' + ion + '/Ion-weighted'
            hfile.require_group(group)
            
           # add attributes
            for(key, value) in projection.items():
                try:
                    hfile[group].attrs[key] = value
                except:
                    continue
                pix_kms = group + '/' + 'pixel_kms'
                self.WriteVariable(projection["pixel_kms"], varname = pix_kms)
            
            # add ion-weighted variables
            for variable in variables:
                self.WriteVariable(projection["Ion-weighted"][ion][variable], 
                                   group+'/' + variable)
            # check for simIon 
            if 'SimIon-weighted' in projection.keys():
                if ion in projection['SimIon-weighted'].keys():
                    group = groupname + '/' + element + '/' + ion + '/SimIon-weighted'
                    hfile.require_group(group)
                    
                # add attributes
                    for(key, value) in projection.items():
                        try:
                            hfile[group].attrs[key] = value
                        except:
                            continue
                        pix_kms = group + '/' + 'pixel_kms'
                        self.WriteVariable(projection["pixel_kms"], varname = pix_kms)
                    
                    # add ion-weighted variables
                    for variable in variables:
                        self.WriteVariable(projection["SimIon-weighted"][ion][variable], 
                                        group+'/' + variable)
                

            
            # optical depth-weighted properties
            group = groupname + '/' + element + '/' + ion + '/optical depth-weighted'
            hfile.require_group(group)
            
           # add attributes
            for(key, value) in projection.items():
                try:
                    hfile[group].attrs[key] = value
                except:
                    continue
                pix_kms = group + '/' + 'pixel_kms'
                self.WriteVariable(projection["pixel_kms"], varname = pix_kms)

            # add transition properties for each transition
            lambda0 = projections["Projection"]["Ion-weighted"][ion]["lambda0"]
            fvalue  = projections["Projection"]["Ion-weighted"][ion]["f-value"]
            #
            info     = {'Vardescription':'Laboratory wavelength', 'CGSConvertsionFactor':1e-8, 'h-scale-exponent':0, 'aexp-scale-exponent':0}
            variable = {'Value': lambda0, 'Info':info}
            self.WriteVariable(variable, varname = group + '/lambda0')
            info     = {'Vardescription':'transition strength', 'CGSConvertsionFactor':1, 'h-scale-exponent':0, 'aexp-scale-exponent':0}
            variable = {'Value': fvalue, 'Info':info}    
            self.WriteVariable(variable, varname = group + '/f-value')            
            # check for simions
            if 'SimIons' in opticaldepth.keys():
                if transition in opticaldepth['SimIons'].keys():
                    for variable in variables:
                        self.WriteVariable(opticaldepth['SimIons'][transition][variable], 
                                group+'/SimIons/' + variable)
                    self.WriteVariable(opticaldepth['SimIons'][transition]['Optical depths'], 
                            group+'/SimIons'+'/Optical depths')

            # add optical dept-weighted variables
            for variable in variables:
                self.WriteVariable(opticaldepth[transition][variable], 
                                   group+'/' + variable)
            self.WriteVariable(opticaldepth[transition]['Optical depths'], 
                                   group+'/Optical depths')
            
             
            
        
        hfile.close()