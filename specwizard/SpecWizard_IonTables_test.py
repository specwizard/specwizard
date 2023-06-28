
import sys
import h5py
import re
import numpy as np 
import scipy.interpolate as interpolate


class IonTables:
    def __init__(self, specparams = None):
        self.specparams   = specparams
        self.iontable_info = specparams['ionparams']
        
    def ReadIonizationTable(self, ion='H I'):
        
        '''
        This will read the ionization tables for the interpolation. The specifics of the table (path,type) are provided by
        initializing the class by providing the output from the Build_Input.
        Input: String of ion to be calculated e.g 'HeII'
        Output: Tuple with the contents of the ionization table for that ion. 
                if iontable_info['specwizard_cloudy'] -> ((redshift,Log_temperature,Log_nH),Log_Abundance)
                if iontable_info['ploeckinger']       -> ((redshift,Log_temperature,Log_(Z/Zsol),Log_nH),Log_Abundance)
        '''     
        iontable_info = self.iontable_info
        if iontable_info["table_type"] == 'specwizard_cloudy':
            # check whether ion exists
            if not ion in iontable_info['ions-available']:
                print("Ion not found")
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
        

    def RomanNumeral(self, i):
        '''
        Roman numerals to decimal. 
        Input: Roman string e.g 'IV'
        Output: Decimal int e.g 4 
        '''
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
        values[values<min_val] = min_val
        values[values>max_val] = max_val
        return values
    
    def IonAbundance(self, redshift=0.0, nH_density=1.0, temperature=1e4, Z=1, ion = 'H I'):
        ''' Return the fractional abundance of a given ion
              
        Compute and return the ionization levels for the given element
        Inputs:
          redshift     : current redshift
          density:     : proper hydrogen density in particles / cm^3 of the gas
          temperature  : temperature of the gas in K
          metallicity  : total metallicity of the particles (only used for ploeckinger tables)
          ionname      : name of ion,e.g. 'H I' for neutral hydrogen
        Output: log_10 n_ion/n_element
        '''
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
            Z          /= Zsol
            # The tables treat zero as log10(1e-50)
            Z[Z==0]    = 1e-50
            TInterpol  = np.log10(temperature)
            Tinterpol  = self.SetLimitRange(TInterpol,table_LogTs.min(),table_LogTs.max())
            nHInterpol = np.log10(nH_density)
            nHInterpol  = self.SetLimitRange(nHInterpol,table_LognHs.min(),table_LognHs.max())
            zInterpol  = redshift
            Zinterpol  = np.log10(Z)
            pts        = np.column_stack((zInterpol, TInterpol, Zinterpol, nHInterpol))
            result     = interpolate.interpn((table_z, table_LogTs, table_LogZs, table_LognHs), table, pts, method='linear', bounds_error=False, fill_value=None)
            return result
        


            return SimulationIonFractions    