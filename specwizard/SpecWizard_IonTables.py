#%load SpecWizard_IonTables.py
#%%writefile SpecWizard_IonTables.py

import sys
import h5py
import re
import roman
import os
import numpy as np 
import scipy.interpolate as interpolate

CHIMES_DICT = {"elec": 0,
                    "HI": 1,
                    "HII": 2,
                    "Hm": 3,
                    "HeI": 4,
                    "HeII": 5,
                    "HeIII": 6,
                    "CI": 7,
                    "CII": 8,
                    "CIII": 9,
                    "CIV": 10,
                    "CV": 11,
                    "CVI": 12,
                    "CVII": 13,
                    "Cm": 14,
                    "NI": 15,
                    "NII": 16,
                    "NIII": 17,
                    "NIV": 18,
                    "NV": 19,
                    "NVI": 20,
                    "NVII": 21,
                    "NVIII": 22,
                    "OI": 23,
                    "OII": 24,
                    "OIII": 25,
                    "OIV": 26,
                    "OV": 27,
                    "OVI": 28,
                    "OVII": 29,
                    "OVIII": 30,
                    "OIX": 31,
                    "Om": 32,
                    "NeI": 33,
                    "NeII": 34,
                    "NeIII": 35,
                    "NeIV": 36,
                    "NeV": 37,
                    "NeVI": 38,
                    "NeVII": 39,
                    "NeVIII": 40,
                    "NeIX": 41,
                    "NeX": 42,
                    "NeXI": 43,
                    "MgI": 44,
                    "MgII": 45,
                    "MgIII": 46,
                    "MgIV": 47,
                    "MgV": 48,
                    "MgVI": 49,
                    "MgVII": 50,
                    "MgVIII": 51,
                    "MgIX": 52,
                    "MgX": 53,
                    "MgXI": 54,
                    "MgXII": 55,
                    "MgXIII": 56,
                    "SiI": 57,
                    "SiII": 58,
                    "SiIII": 59,
                    "SiIV": 60,
                    "SiV": 61,
                    "SiVI": 62,
                    "SiVII": 63,
                    "SiVIII": 64,
                    "SiIX": 65,
                    "SiX": 66,
                    "SiXI": 67,
                    "SiXII": 68,
                    "SiXIII": 69,
                    "SiXIV": 70,
                    "SiXV": 71,
                    "SI": 72,
                    "SII": 73,
                    "SIII": 74,
                    "SIV": 75,
                    "SV": 76,
                    "SVI": 77,
                    "SVII": 78,
                    "SVIII": 79,
                    "SIX": 80,
                    "SX": 81,
                    "SXI": 82,
                    "SXII": 83,
                    "SXIII": 84,
                    "SXIV": 85,
                    "SXV": 86,
                    "SXVI": 87,
                    "SXVII": 88,
                    "CaI": 89,
                    "CaII": 90,
                    "CaIII": 91,
                    "CaIV": 92,
                    "CaV": 93,
                    "CaVI": 94,
                    "CaVII": 95,
                    "CaVIII": 96,
                    "CaIX": 97,
                    "CaX": 98,
                    "CaXI": 99,
                    "CaXII": 100,
                    "CaXIII": 101,
                    "CaXIV": 102,
                    "CaXV": 103,
                    "CaXVI": 104,
                    "CaXVII": 105,
                    "CaXVIII": 106,
                    "CaXIX": 107,
                    "CaXX": 108,
                    "CaXXI": 109,
                    "FeI": 110,
                    "FeII": 111,
                    "FeIII": 112,
                    "FeIV": 113,
                    "FeV": 114,
                    "FeVI": 115,
                    "FeVII": 116,
                    "FeVIII": 117,
                    "FeIX": 118,
                    "FeX": 119,
                    "FeXI": 120,
                    "FeXII": 121,
                    "FeXIII": 122,
                    "FeXIV": 123,
                    "FeXV": 124,
                    "FeXVI": 125,
                    "FeXVII": 126,
                    "FeXVIII": 127,
                    "FeXIX": 128,
                    "FeXX": 129,
                    "FeXXI": 130,
                    "FeXXII": 131,
                    "FeXXIII": 132,
                    "FeXXIV": 133,
                    "FeXXV": 134,
                    "FeXXVI": 135,
                    "FeXXVII": 136,
                    "H2": 137,
                    "H2p": 138,
                    "H3p": 139,
                    "OH": 140,
                    "H2O": 141,
                    "C2": 142,
                    "O2": 143,
                    "HCOp": 144,
                    "CH": 145,
                    "CH2": 146,
                    "CH3p": 147,
                    "CO": 148,
                    "CHp": 149,
                    "CH2p": 150,
                    "OHp": 151,
                    "H2Op": 152,
                    "H3Op": 153,
                    "COp": 154,
                    "HOCp": 155,
                    "O2p": 156}
class IonTables:
    """
    Reads and interpolates ionization tables for calculating ion fractions.

    Attributes:
        specparams (dict): Wizard dictionary produced by the `BuildInput` routine.
        iontable_info (dict): Metadata and configuration for the ionization tables.
    """

    def __init__(self, specparams = None):
        """
        Initializes the IonTables class.

        Args:
            specparams (dict): Dictionary containing parameters, including the ionization table configuration.
        """        
        self.specparams   = specparams
        self.iontable_info = specparams['ionparams']
       
        
    def ReadIonizationTable(self, ion='H I'):
        
        """
        Reads the ionization table for the specified ion.

        Args:
            ion (str): Name of the ion to be calculated, e.g., 'HeII'.

        Returns:
            tuple: Ionization table contents.
                - If `table_type` is 'specwizard_cloudy': ((redshift, Log_temperature, Log_nH), Log_Abundance).
                - If `table_type` is 'ploeckinger': ((redshift, Log_temperature, Log_(Z/Zsol), Log_nH), Log_Abundance).

        Raises:
            FileNotFoundError: If the required ionization table file does not exist.
            ValueError: If the ion is not available in the table configuration.
        """
        iontable_info = self.iontable_info
        if iontable_info["table_type"] == 'specwizard_cloudy':
            # check whether ion exists
            ion_available_names = np.array([ion for _,ion in iontable_info['ions-available']])
            if not ion in ion_available_names:
                print("Ion",ion," is not found")
                sys.exit(-1)
            fname = iontable_info['iondir'] + '/' + ion + '.hdf5'
            hf    = h5py.File(fname,  "r")
                
            # cloudy parameters
            UVB   = hf.attrs.get("UVB")
            compo = hf.attrs.get("Composition")
            
            # tables
            z            = hf['z'][...]
            LogT         = hf['LogT'][...]
            LognH        = hf['LognH'][...]
            LogAbundance = hf['LogAbundance'][...]
            #
            hf.close()
            #
            return ((z, LogT, LognH), LogAbundance)
        
        if iontable_info["table_type"] == 'ploeckinger':

            file       = iontable_info["iondir"] + '/' + iontable_info["fname"]
            hf         = h5py.File(file,  "r")
            # read temperatre bins
            LogT       = hf["TableBins/TemperatureBins"][...]
            LognH      = hf["TableBins/DensityBins"][...]
            z          = hf["TableBins/RedshiftBins"][...]
            LogZ       = hf["TableBins/MetallicityBins"][...]
            # read and format the correct element
            # 
            element    = (''.join(re.split('I|  |V|X|', ion))).strip()
            
            #print("element to read = {}".format(element))
            # We create a dictionary between the element names (as provided to from the user) and
            # how they are called in the Ploeckinger tables e,g H -> 01hydrogen. 
            srtnames   = np.array(hf['ElementNamesShort'][...],dtype=str)
            ionnames   = np.array(list(hf['Tdep/IonFractions'].keys()),dtype=str)
            tablenames = dict(zip(srtnames,ionnames))            
            tablename  = tablenames[element]
            
 #           print("table name = ", tablename)
            # We use self.RomanNumeral to create a dictionary between Roman and decimal numerals
            # So we can know that the user input O VI -> O 6  
            dec_and_rom = np.array([[i,self.RomanNumeral(i)] for i in range(30)])
            rom2dec_dict = dict(zip(dec_and_rom[:,1],dec_and_rom[:,0]))
            rom_num       = ((ion).replace(element,"")).strip()
            ionlevel     = int(rom2dec_dict[rom_num]) - 1
 #           print("ion = ", ion, 'level =', ionlevel)
            #
            LogAbundance = hf['Tdep/IonFractions/' + tablename][:,:,:,:,ionlevel]
 #           LogAbundance = hf['Tdep/HydrogenFractionsVol'][:,:,:,:,0]
            #
            hf.close()
            #
            
            return ((z, LogT, LognH, LogZ), LogAbundance)



        if iontable_info["table_type"] == 'chimes':
            #CHIMES needs romal numerals to be removed and spaces to be removed to match the keys in the CHIMES_DICT
            ion_key = ion.replace(' ', '')
            

            element_index = CHIMES_DICT[ion_key]
            iondir = iontable_info["iondir"]

            # Build a redshift-resolved table by stacking all z*.hdf5 files in iondir.
            redshift_pattern = re.compile(r"^z(?P<z>\d+(?:\.\d+)?)(?:_eqm)?\.hdf5$")
            redshift_files = []
            for onefile in os.listdir(iondir):
                match = redshift_pattern.match(onefile)
                if match:
                    redshift_files.append((float(match.group("z")), os.path.join(iondir, onefile)))

            if len(redshift_files) == 0:
                raise FileNotFoundError(
                    f"No CHIMES redshift files found in {iondir}. Expected files like z0.000.hdf5 or z0.000_eqm.hdf5"
                )

            redshift_files.sort(key=lambda x: x[0])
            z = np.array([item[0] for item in redshift_files])

            first_file = redshift_files[0][1]
            with h5py.File(first_file, "r") as hf0:
                LogT = hf0["TableBins/Temperatures"][:]
                LognH = hf0["TableBins/Densities"][:]
                LogZ = hf0["TableBins/Metallicities"][:]

            # Final table ordering matches ploeckinger interpolation usage: (z, T, Z, nH).
            LogAbundance = np.empty((len(z), len(LogT), len(LogZ), len(LognH)))

            for i, (_, filepath) in enumerate(redshift_files):
                with h5py.File(filepath, "r") as hf:
                    file_LogT = hf["TableBins/Temperatures"][:]
                    file_LognH = hf["TableBins/Densities"][:]
                    file_LogZ = hf["TableBins/Metallicities"][:]                  
                    
                     #check for every file that the bins are the same, 
                     #otherwise interpolation will be wrong and we will have to interpolate to a common grid firsties"][:]
                    if not (np.array_equal(file_LogT, LogT) and np.array_equal(file_LognH, LognH) and np.array_equal(file_LogZ, LogZ)):
                        raise ValueError(f"Inconsistent CHIMES table bins in {filepath}")

                    # File layout is (T, nH, Z, species). Move to ( to match original ordering in "ploeckinger table", Z, nH).
                    LogAbundance[i] = np.transpose(hf["Abundances"][:, :, :, element_index], (0, 2, 1))

            result = ((z, LogT, LognH, LogZ), LogAbundance)
            
            return result
        
        if iontable_info["table_type"] == 'cloudy_hm01':
            ion_element_short = (''.join(re.split('I|  |V|X|', ion))).strip().lower()
            ion_number = roman.fromRoman(ion[len((''.join(re.split('I|  |V|X|', ion))).strip()):].strip())
            ion_file = f"{ion_element_short}{ion_number}.hdf5"
            # check whether ion exists
#            if f"{ion_element_short}{ion_number}" not in iontable_info['ions-available']:
#                print(f"Ion {ion_element_short}{ion_number} not found")
#                sys.exit(-1)
            fname = os.path.join(iontable_info['iondir'], ion_file)
            hf    = h5py.File(fname,  "r")

        #            # cloudy parameters
        #            UVB   = hf.attrs.get("UVB")
        #            compo = hf.attrs.get("Composition")

            # tables
            z            = hf['redshift'][...]
            LogT         = hf['logt'][...]
            LognH        = hf['logd'][...]
            LogAbundance = np.log10(hf['ionbal'][...])
            #
            hf.close()
            #
            return ((z, LogT, LognH), LogAbundance)


        
    def RomanNumeral(self, i):
        """
        Converts a decimal number to a Roman numeral.

        Args:
            i (int): Decimal number.

        Returns:
            str: Roman numeral representation.
        """
        result = ''
        # number of Xs
        nX = int(i/10)
        for c in np.arange(nX):
            result += 'X'
        j  = i%10
        if j < 4:
            for c in np.arange(j):
                result += 'I'
        elif j==4:
            result += 'IV'
        elif j == 5:
            result += 'V'
        elif j == 9:
            result += 'IX'
        else:
            k = j - 5
            result += 'V'
            for c in np.arange(k):
                result += 'I'
        return result
    
    def SetLimitRange(self,values,min_val,max_val):
        """
        Clamps values to a specified range.

        Args:
            values (np.ndarray): Array of values to limit.
            min_val (float): Minimum allowed value.
            max_val (float): Maximum allowed value.

        Returns:
            np.ndarray: Array with values limited to the specified range.
        """
        values[values<min_val] = min_val
        values[values>max_val] = max_val
        return values
    
    def IonAbundance(self, redshift=None, nH_density=None, temperature=None, metal_fraction=None, ion = None):
        """
        Computes the fractional abundance of a specified ion.

        Args:
            redshift (float): Redshift of the particles.
            nH_density (float): Proper hydrogen density (particles/cm³).
            temperature (float): Gas temperature in Kelvin.
            metal_fraction (float): Total metallicity of the particles (for 'ploeckinger' tables).
            ion (str): Name of the ion (e.g., 'H I').

        Returns:
            float: Ion fraction (log10 n_ion / n_element).

        Raises:
            ValueError: If required data is missing or interpolation fails.
        """
        iontable_info = self.iontable_info

        if iontable_info["table_type"] == 'specwizard_cloudy':
            # read the table
            (table_z, table_LogTs, table_LognHs), table = self.ReadIonizationTable(ion=ion)
            #
            TInterpol  = np.log10(temperature)
            Tinterpol  = self.SetLimitRange(TInterpol,table_LogTs.min(),table_LogTs.max())
            nHInterpol = np.log10(nH_density)
            nHInterpol  = self.SetLimitRange(nHInterpol,table_LognHs.min(),table_LognHs.max())
            zInterpol  = redshift
            pts        = np.column_stack((zInterpol, TInterpol, nHInterpol))
            result     = interpolate.interpn((table_z, table_LogTs, table_LognHs), table, pts, method='linear', bounds_error=False, fill_value=None)
            return result
        if iontable_info["table_type"] == 'ploeckinger':
            ((table_z, table_LogTs, table_LognHs, table_LogZs), table) = self.ReadIonizationTable(ion=ion)
            # We divide by the Solar metallicity since the tables are in log10(Z/Zsol)
            Zsol       = 0.013371374
            Z          = metal_fraction/Zsol
            # The tables treat zero as log10(1e-50)
            TInterpol  = np.log10(temperature)
            Tinterpol  = self.SetLimitRange(TInterpol,table_LogTs.min(),table_LogTs.max())
            nHInterpol = np.log10(nH_density)
            nHInterpol  = self.SetLimitRange(nHInterpol,table_LognHs.min(),table_LognHs.max())
            zInterpol  = redshift
            Zinterpol  = np.log10(Z)
#            Zinterpol[Zinterpol==-34.00] = -50.00
            Zinterpol  = self.SetLimitRange(Zinterpol,table_LogZs.min(),table_LogZs.max())
            pts        = np.column_stack((zInterpol, TInterpol, Zinterpol, nHInterpol))
            result     = interpolate.interpn((table_z, table_LogTs, table_LogZs, table_LognHs), table, pts, method='linear', bounds_error=False, fill_value=None)
            return result

        if iontable_info["table_type"] == 'chimes':
            ((table_z, table_LogTs, table_LognHs, table_LogZs), table) = self.ReadIonizationTable(ion=ion)
            # We divide by the Solar metallicity since the tables are in log10(Z/Zsol)
            Zsol       = 0.013371374
            Z          = metal_fraction/Zsol
            # The tables treat zero as log10(1e-50)
            TInterpol  = np.log10(temperature)
            Tinterpol  = self.SetLimitRange(TInterpol,table_LogTs.min(),table_LogTs.max())
            nHInterpol = np.log10(nH_density)
            nHInterpol  = self.SetLimitRange(nHInterpol,table_LognHs.min(),table_LognHs.max())
            zInterpol  = redshift
            Zinterpol  = np.log10(Z)
#            Zinterpol[Zinterpol==-34.00] = -50.00
            Zinterpol  = self.SetLimitRange(Zinterpol,table_LogZs.min(),table_LogZs.max())
            pts        = np.column_stack((zInterpol, TInterpol, Zinterpol, nHInterpol))
            result     = interpolate.interpn((table_z, table_LogTs, table_LogZs, table_LognHs), table, pts, method='linear', bounds_error=False, fill_value=None)
            return result
        
        if iontable_info["table_type"] == 'cloudy_hm01':
            # Read the ionization table
            ((table_redshifts, table_log_temps, table_log_hydrogen_number_density), table_abundances) = self.ReadIonizationTable(ion=ion)

#            # Convert metalicities
#            Z /= 0.013371374 # We divide by the Solar metallicity since the tables are in log10(Z/Zsol)
#            Z[Z==0] = 1e-50 # The tables treat zero as log10(1e-50)#TODO: is this true for the older tables???

            # Set the limits for the 
            TInterpol  = self.SetLimitRange(np.log10(temperature),
                                            table_log_temps.min(), table_log_temps.max())
            nHInterpol = self.SetLimitRange(np.log10(nH_density),
                                            table_log_hydrogen_number_density.min(), table_log_hydrogen_number_density.max())
#            Zinterpol = np.log10(np.where(Z > 10e-35, Z, 10e-35))#TODO: handle 0 metalicities properly

            # Define where to sample the 4D grid at
            pts = np.column_stack((nHInterpol, TInterpol, redshift))
            
            # Return the samples
            return interpolate.interpn((table_log_hydrogen_number_density, table_log_temps, table_redshifts), table_abundances,
                                       pts,
                                       method='linear', bounds_error=False, fill_value=None)
        
        
        
        

            return SimulationIonFractions    
