import numpy as np
from os.path import exists
from os import remove
import copy
import h5py
import unyt
#
#from SpecWizard_Elements import Elements
#from SpecWizard_Input import ReadSnap
#from SpecWizard_SplineInterpolation import ColumnTable
#from SpecWizard_IonizationBalance import IonizationBalance
#from SpecWizard_ProjectData import SightLineProjection

    
class OpticalDepth_IO:
    ''' Methods to write and read files containing spectra and optical depth '''
    def __init__(self, wizard=None, create=True):
        self.dirname    = wizard["Output"]["directory"]
        self.fname      = wizard["Output"]["fname"]
        self.wizard = wizard
        self.nspec      = 0
        try:
            self.header     = wizard["Header"]
        except:
            print("Warning: No header information found in wizard. Output file will be missing header information.")
            self.header = None
        
       
        
        # if create: create file and add header. This is the write mode
        if create:
            file = self.dirname + self.fname
            self.ions       = wizard["ionparams"]["Ions"]
            element_names   = np.array([ self.ions[i][0] for i in range(len(self.ions))])
            self.elements   = element_names
            if exists(file):
                print("Removing file ", file)
                remove(file)
            # write header only if it exists
            if self.header is not None:
                self.WriteHeader()
        #if not create: this is the read mode
        
            
    def ReadGroup(self, groupname='Header'):
        # read all entries for this particular hdf5 group
        hfile = h5py.File(self.dirname + self.fname, "r")
        group = hfile[groupname]
        grp_dict = {}
        for k in group.attrs.keys():
            grp_dict[k] = group.attrs[k]

        hfile.close()
        return dict(grp_dict)

    def ReadVariable(self, varname='LOS1/Density'):
        # read one dataset and all of its attributes
        hfile = h5py.File(self.dirname + self.fname, "r")
        dset = hfile[varname]

        values = dset[...]
        info = {}
        for k in dset.attrs.keys():
            info[k] = dset.attrs[k]

        hfile.close()
        return {'Value': values, 'Info': dict(info)}

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
        
        # create and write dataset
        Values = variable["Value"]
        # Unyt does not like to read unyt_quantities
        if isinstance(Values, unyt.unyt_quantity):
            Values = unyt.unyt_array(Values)

        #hfile.create_dataset(varname, data=Values)
        just_varname = varname.split('/')[-1]
        dir_to_group = varname.replace(just_varname,'')
        Values.write_hdf5(self.dirname + self.fname,dataset_name=just_varname, group_name=dir_to_group)
        hfile  = h5py.File(self.dirname + self.fname, "a")
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

    def ReadHeader(self):
        """
        Read header attributes and fields
        """
        header = {}

        try:
            header['attributes'] = self.ReadGroup("Header")
        except:
            header['attributes'] = {}

        header['datasets'] = {}

        # discover datasets inside Header
        with h5py.File(self.dirname + self.fname, "r") as f:
            if "Header" in f:
                for field in f["Header"]:
                    try:
                        header['datasets'][field] = self.ReadVariable(f"Header/{field}")
                    except:
                        continue

        return header


    def write_shortspectra_to_file(self, projections):
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
        od_variables = variables + ['TotalIonColumnDensity', 'HydrogenDensities', 'Metallicities']
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
            group = groupname + '/' + element + '/' + ion + '/Ion optical depth-weighted'
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
            lambda0 = unyt.unyt_array(projections["Projection"]["Ion-weighted"][ion]["lambda0"])
            fvalue  = unyt.unyt_array(projections["Projection"]["Ion-weighted"][ion]["f-value"])
            #
            info     = {'Vardescription':'Laboratory wavelength', 'h-scale-exponent':0, 'aexp-scale-exponent':0}
            variable = {'Value': lambda0, 'Info':info}
            self.WriteVariable(variable, varname = group + '/lambda0')
            info     = {'Vardescription':'transition strength', 'h-scale-exponent':0, 'aexp-scale-exponent':0}
            variable = {'Value': fvalue, 'Info':info}    
            self.WriteVariable(variable, varname = group + '/f-value') 

            # add optical dept-weighted variables
            for variable in od_variables:
                if variable in opticaldepth[transition].keys():
                    self.WriteVariable(opticaldepth[transition][variable], 
                                       group+'/' + variable)
            self.WriteVariable(opticaldepth[transition]['Optical depths'], 
                                   group+'/Optical depths')  

            # check for simions, if so please add simion optical depth-weighted variables as well
            if 'SimIons' in opticaldepth.keys():
                
                if transition in opticaldepth['SimIons'].keys():
                    group = groupname + '/' + element + '/' + ion + '/SimIon optical depth-weighted'
                    hfile.require_group(group)
                    for variable in od_variables:
                        if variable in opticaldepth['SimIons'][transition].keys():
                            self.WriteVariable(opticaldepth['SimIons'][transition][variable], 
                                    group + '/' + variable)
                    self.WriteVariable(opticaldepth['SimIons'][transition]['Optical depths'], 
                            group+'/Optical depths')
                    # add transition properties for each transition, copy from the ion weighted properties
                    lambda0 = unyt.unyt_array(projections["Projection"]["Ion-weighted"][ion]["lambda0"])
                    fvalue  = unyt.unyt_array(projections["Projection"]["Ion-weighted"][ion]["f-value"])
                    #
                    info     = {'Vardescription':'Laboratory wavelength', 'h-scale-exponent':0, 'aexp-scale-exponent':0}
                    variable = {'Value': lambda0, 'Info':info}
                    self.WriteVariable(variable, varname = group + '/lambda0')
                    info     = {'Vardescription':'transition strength', 'h-scale-exponent':0, 'aexp-scale-exponent':0}
                    variable = {'Value': fvalue, 'Info':info}    
                    self.WriteVariable(variable, varname = group + '/f-value')

                    # add attributes, also copy from the ion weighted properties
                    for(key, value) in projection.items():
                        try:
                            hfile[group].attrs[key] = value
                        except:
                            continue
                        pix_kms = group + '/' + 'pixel_kms'
                        self.WriteVariable(projection["pixel_kms"], varname = pix_kms)
                

            
            
             
            
        
        hfile.close()
    
    def read_shortspectra_from_file(self):
        """
        Clean short spectra reader with structure:

        result["Header"]
        result["Data"][element][ion][subgroup][field]
        """

        import h5py

        file = self.dirname + self.fname

        result = {
            "Header": self.ReadHeader(),
            "Data": {}
        }

        with h5py.File(file, "r") as f:

            for los in f.keys():

                if not los.startswith("LOS_"):
                    continue

                data_dict = {}

                # -------------------------
                # optional: store LOS attributes at top level
                # -------------------------
                data_dict["attributes"] = self.ReadGroup(los)

                if f"{los}/Box_kms" in f:
                    data_dict["Box_kms"] = self.ReadVariable(f"{los}/Box_kms")

                # -------------------------
                # element loop
                # -------------------------
                for element in f[los].keys():

                    if element == "Box_kms":
                        continue

                    element_path = f"{los}/{element}"
                    element_dict = {}

                    # -------------------------
                    # loop over element content
                    # -------------------------
                    for item in f[element_path].keys():

                        item_path = f"{element_path}/{item}"

                        #element weighted properties
                        if item == "Element-weighted":

                            ew_dict = {}
                            for field in f[item_path].keys():
                                full_path = f"{item_path}/{field}"
                                ew_dict[field] = self.ReadVariable(full_path)

                            element_dict["Element-weighted"] = ew_dict
                            continue

                        # ion properties
                        ion = item
                        ion_path = item_path
                        ion_dict = {}

                        for subgroup in f[ion_path].keys():

                            subgroup_path = f"{ion_path}/{subgroup}"
                            subgroup_obj = f[subgroup_path]

                            sub_dict = {}

                            for field in subgroup_obj.keys():

                                full_path = f"{subgroup_path}/{field}"
                                obj = f[full_path]
                                sub_dict[field] = self.ReadVariable(full_path)


                            ion_dict[subgroup] = sub_dict

                        element_dict[ion] = ion_dict

                    data_dict[element] = element_dict

                
                result["Data"][los] = data_dict

                # assume one LOS per file
                break

        return result



    def write_longspectra_to_file(self, projections):
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
        od_variables = variables + ['TotalIonColumnDensity', 'HydrogenDensities', 'Metallicities']
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
                    for variable in od_variables:
                        if variable in opticaldepth['SimIons'][transition].keys():
                            self.WriteVariable(opticaldepth['SimIons'][transition][variable], 
                                    group+'/SimIons/' + variable)
                    self.WriteVariable(opticaldepth['SimIons'][transition]['Optical depths'], 
                            group+'/SimIons'+'/Optical depths')

            # add optical dept-weighted variables
            for variable in od_variables:
                if variable in opticaldepth[transition].keys():
                    self.WriteVariable(opticaldepth[transition][variable], 
                                       group+'/' + variable)
            self.WriteVariable(opticaldepth[transition]['Optical depths'], 
                                   group+'/Optical depths')
            
             
            
        
        hfile.close()

    def write_fullspectrum_to_file(self, long_spectra,header_dict=None):
        """
        Save the full long spectrum (from do_long_spectra) to a HDF5 file.
        """

        hfile = h5py.File(self.dirname + self.fname, "a")

        base_group = "FullSpectrum"
        hfile.require_group(base_group)

        # --- Write global grids ---
        self.WriteVariable(
            {'Value': long_spectra['velocities'],
            'Info': {'Vardescription': 'velocity grid'}},
            f'{base_group}/Velocities'
        )

        self.WriteVariable(
            {'Value': long_spectra['wavelengths'],
            'Info': {'Vardescription': 'wavelength grid'}},
            f'{base_group}/Wavelengths'
        )

        # optional header information for the full spectrum, not yet implemented
        if header_dict is not None:
            header_group = f"{base_group}/Header"
            hfile.require_group(header_group)
            for key, value in header_dict.items():
                hfile[header_group].attrs[key] = value

        
        ions_group = f"{base_group}"
        hfile.require_group(ions_group)

        for transition in long_spectra['Ions']:
            element, ion = transition

            element_group = f"{ions_group}/{element}"
            hfile.require_group(element_group)

            ion_group = f"{element_group}/{ion}"
            hfile.require_group(ion_group)

            ion_data = long_spectra['Ions'][transition]

            # --- transition properties ---
            self.WriteVariable(
                {'Value': ion_data['lambda0'],
                'Info': {'Vardescription': 'Laboratory wavelength',
                        'h-scale-exponent': 0,
                        'aexp-scale-exponent': 0}},
                f"{ion_group}/lambda0"
            )

            self.WriteVariable(
                {'Value': ion_data['f-value'],
                'Info': {'Vardescription': 'transition strength',
                        'h-scale-exponent': 0,
                        'aexp-scale-exponent': 0}},
                f"{ion_group}/f-value"
            )

            # --- other fields ---
            for field, data in ion_data.items():

                if field in ['lambda0', 'f-value']:
                    continue

                dataset_name = field.replace(" ", "_")
                full_path = f"{ion_group}/{dataset_name}"

                if isinstance(data, dict) and 'Value' in data:
                    self.WriteVariable(data, full_path)
                else:
                    self.WriteVariable(
                        {'Value': data,
                        'Info': {'Vardescription': field}},
                        full_path
                    )

        hfile.close()

    def read_fullspectrum_from_file(self):

        result = {
            "Header": self.ReadHeader(),
            "Data": {}
        }

        file = self.dirname + self.fname
        base = "FullSpectrum"

        with h5py.File(file, "r") as f:

            # global values
            result["Data"]["Velocities"] = self.ReadVariable(f"{base}/Velocities")
            result["Data"]["Wavelengths"] = self.ReadVariable(f"{base}/Wavelengths")

            # element loop
            for element in f[base].keys():

                if element in ["Velocities", "Wavelengths", "Header"]:
                    continue

                result["Data"][element] = {}
                element_path = f"{base}/{element}"

                # Ion loop
                for ion in f[element_path].keys():

                    ion_path = f"{element_path}/{ion}"
                    ion_dict = {}

                    
                    for field in f[ion_path].keys():

                        try:
                            ion_dict[field] = self.ReadVariable(f"{ion_path}/{field}")
                        except:
                            continue

                    result["Data"][element][ion] = ion_dict

        return result