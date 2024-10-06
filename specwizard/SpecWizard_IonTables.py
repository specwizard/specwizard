#%load SpecWizard_IonTables.py
#%%writefile SpecWizard_IonTables.py

import sys
import h5py
import re
import roman
import os
import numpy as np 
import scipy.interpolate as interpolate


class IonTables:
    """
    This class reads and interpolates the ionization tables that will be used for the calculation of ion fractions. 

    Parameters
    ----------

    specparams : dictionary
        The Wizard dictionary produced by the BuildInput routine. 
    """

    def __init__(self, specparams = None):
        self.specparams   = specparams
        self.iontable_info = specparams['ionparams']
        
    def ReadIonizationTable(self, ion='H I'):
        
        """
        This will read the ionization tables for the interpolation. The specifics of the table (path,type) are provided by
        initializing the class by providing the output from the Build_Input.
	
	Parameters
	----------
        
        ion : str
	      Name of ion to be calculated e.g 'HeII'

	Returns
	-------
        
	out : tuple
		tuple with the contents of the ionization table for that ion.
                if 'specwizard_cloudy' -> ((redshift,Log_temperature,Log_nH),Log_Abundance)
                if 'ploeckinger'       -> ((redshift,Log_temperature,Log_(Z/Zsol),Log_nH),Log_Abundance)
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
        Converts decimal to roman numerals . 

        Parameters
	----------

        i : int
	Decimal number

	Returns
	-------
	
        out : str
        Roman string e.g 'IV'
         
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
        Set values < min_val to min_val and values > min_max to max_val. This to avoid extrapolation. 

        Parameters
        ----------

        values : np.array(float)
            Array with the values we want to limit

        min_val : float 
            Minimal value range set by the table. 

        max_val : float 
            Maximal value range set by the table. 

        Returns 
        -------
        out : np.array(float)
            array with values bounded by the table.  
	"""
        values[values<min_val] = min_val
        values[values>max_val] = max_val
        return values
    
    def IonAbundance(self, redshift=None, nH_density=None, temperature=None, metal_fraction=None, ion = None):
        """ 
	Return the fractional abundance of a given ion. This done by interpolation of the ionization tables with particle data from simulations. 
 	             
        Parameters
	----------

        redshift     : float
        Redshift of the particles. 

        density:     : float
	Proper hydrogen density of gas  particles in $\mathrm{cm}^{-3}$

        temperature  : float 
	Temperature of  gas particles in $\mathrm{K}$

        metal_fraction  : float
	Total metallicity of the particles (only used for ploeckinger tables)
          
	ionname      : str
	Name of ion,e.g. 'H I' for neutral hydrogen

	Returns
	-------

        out : float
	Ion fraction $\log_{10} \mathrm{n_{ion}/n_{element}}$
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
