import numpy as np
import requests
from bs4 import BeautifulSoup
import h5py
import traceback
import mendeleev
from mendeleev.fetch import fetch_table
import roman

class Atomfile:
    """
    This class creates the atom file containing atomic transition data used to generate spectral lines for ions.
    The data is retrieved from the NIST Atomic Spectra Database (ASD): https://physics.nist.gov/PhysRefData/ASD/lines_form.html
    """


    def __init__(self,do_all_elements=False,elements_to_do=None):
        """
        Initializes the Atomfile class by setting up URLs for the NIST database and generating an atomic directory.
        
        Args:
            do_all_elements (bool): If True, includes all elements in the periodic table. If False, includes a default subset.
        """
        self.nist_url = 'https://physics.nist.gov/cgi-bin/ASD/lines1.pl?'

        self.nist_grnd_url = 'https://physics.nist.gov/cgi-bin/ASD/ie.pl?'

    # we generate an atom direcotry with information about the element and ions

        romnumerals    = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 
                    'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI', 'XVII',
                'XVIII', 'XIX', 'XX', 'XXI', 'XXII', 'XXIII', 'XXIV', 'XXV', 'XXVI']
        

        if do_all_elements:
                ptable = fetch_table("elements")
                element_names = list(ptable.name)
                element_short = [ mendeleev.element(name).symbol for name in element_names]
                element_num_states = [ getattr(mendeleev,elmt).atomic_number for elmt in element_short]
		
        elif isinstance(elements_to_do, list):
                element_short  =  elements_to_do
                element_names  =  [ getattr(mendeleev,elmt).name for elmt in element_short]
                element_num_states = [ getattr(mendeleev,elmt).atomic_number for elmt in element_short]

        else:
                element_short  = ['H', 'He', 'C', 'N', 'O', 'Ne', 'Mg', 'Si', 'S', 'Ca', 'Fe']
                element_names =  [ getattr(mendeleev,elmt).name for elmt in element_short]
                element_num_states = [ getattr(mendeleev,elmt).atomic_number for elmt in element_short]
	
        # We add Deuterium 
        element_short.insert(1,'D')
        element_names.insert(1,'Deuterium')
        element_num_states.insert(1,1)
        atom_dictionary = {}

        for i,element in enumerate(element_names):
            atom_dictionary[element] = {}
            atom_dictionary[element]["lines"] = {}
            nstate = element_num_states[i]
            atom_dictionary[element]["symbol"] = element_short[i]
            atom_dictionary[element]["lines"] = []
            atom_dictionary[element]["nstates"] = nstate
            for j in range(nstate):
                atom_dictionary[element]["lines"].append(element_short[i]+" "+roman.toRoman(j+1))

        self.atom_dictionary = atom_dictionary 

    def format_line(self,line):
        """
        Adds a zero to any column in the line that has an empty space.
        
        Args:
            line (list): A list of strings representing a line of data.
            
        Returns:
            list: The formatted line with empty spaces replaced by '0'.
        """
        for i,ln in enumerate(line):
            if len(ln.strip()) == 0:
                line[i] = '0'
            else:
                line[i] = ln.strip()
        return line    
    
    def clean_rel_ints(self,value):
        """
        Cleans intensity values from the NIST database by removing unwanted characters.
        
        Args:
            value (str): Intensity value from the NIST database.
            
        Returns:
            str: Cleaned intensity value as a string.
        """
        rel_ints_markers = ['*',':','-',',','a','bl','b','B','c','d','D','E','N','f','F','e','g','G','H','hfs','h','i','j','l','m','o','O','p','q','r','s','S','I','t','u','w','x','(',')','IV/2',' ','V/2','?']
        for marker in rel_ints_markers:
            value = value.replace(marker,"")
        
        if value == '':
            value = '0'
        return value    
    
    def format_lines_into_array(self,lines):
        """
        Formats query results from the NIST database into a numpy array.
        
        Args:
            lines (list): List of strings representing the lines from the query response.
            
        Returns:
            numpy.ndarray: Formatted array containing the data.
        """
        df = []
        # The first 6 lines are the header 
        for line in lines[6:-1]:
            # Every 5 lines a  line with 98 white spaces
            if len(line.strip(' '))==112:
                continue
            else:
                bleh      = line.split('|')
                # We introduce 0 into values in the table that are empty
                bleh_frmt = self.format_line(bleh)
                # we ignore the last two values since are not relevant. 
                df.append(bleh_frmt[:11])
        llines = np.vstack(df)
        return llines    
    
    
    def check_for_query_errors(self,response):
        """
        Checks the query response for potential errors.
        
        Args:
            response (requests.Response): Response object from the NIST database query.
            
        Raises:
            Exception: If the response indicates no lines are available or an unrecognized ion.
        """
        if 'No lines are available in ASD with the parameters selected' in response.text:
            
            raise Exception("Line not in the NIST database with current settings.")
            
        if 'Unrecognized token.' in response.text:
            
            raise Exception("Unrecognized ion")

    def mask_ground_states(self, ground_state,lines):
        """
        Filters the lines to include only those corresponding to the ground state.
        
        Args:
            ground_state (str): Ground level configuration.
            lines (numpy.ndarray): Array of lines to filter.
            
        Returns:
            numpy.ndarray: Filtered lines corresponding to the ground state.
            
        Raises:
            Exception: If no ground state lines are found.
        """
    
        ground_mask =lines[:,5]==ground_state

        if any(ground_mask):
            return lines[ground_mask]
        else:
            raise Exception("Ground state not found!")



    def get_groundstate(self,ion_to_q):
        """
        Queries the NIST database to get the ground level configuration for a specific ion.
        
        Args:
            ion_to_q (str): Ion to query, formatted as "Element Symbol RomanNumeral" (e.g., "H I", "Fe XI").
            
        Returns:
            str: Ground level configuration.
        """


        url_grnd = self.nist_grnd_url+'spectra='+ion_to_q+'&units=1&format=1&remove_js=on&order=0&conf_out=on&e_out=0&submit=Retrieve+Data' 
        respond = requests.get(url=url_grnd)
        table_post_pre = BeautifulSoup(respond.text, features='html5lib').find('pre')
        lines = table_post_pre.text.strip().split('\n')
        ground_state = lines[3].replace('|','').split()[0]
        
        return ground_state        

    def create_hdf5_from_nist(self,file_name="atom_file.hdf5", wavelength_low_lim=200.0,wavelength_upper_lim=8000.0,wavelength_resolution=0.1,verbose=False):
        """
        Creates an HDF5 file containing atomic transition data queried from the NIST database.
        
        Args:
            file_name (str): Name of the HDF5 file.
            wavelength_low_lim (float): Lower wavelength limit in Ångströms.
            wavelength_upper_lim (float): Upper wavelength limit in Ångströms.
            wavelength_resolution (float): Spectral resolution in Ångströms.
        """
        atom_dictionary = self.atom_dictionary
        self.verbose  = verbose
        self.wavelength_low = wavelength_low_lim
        self.wavelength_high = wavelength_upper_lim
        self.delta_lambda    = wavelength_resolution
        atomh5 = h5py.File(file_name,"a")

        for element_name in atom_dictionary.keys():
            for ion_name in atom_dictionary[element_name]["lines"]:
                if self.verbose:
                    print("Writting "+element_name+ " "+ ion_name+ "...")
                try:
                    lines, line_url= self.get_nist_line(ion2do=ion_name)
                    
                except Exception as exc:

                    #print(traceback.format_exc())
                    if self.verbose:
                        print (ion_name+':', exc)
        
                    continue
        
                lambda0 = lines[:,0].astype(float)
                fvalue  = lines[:,3].astype(float)
                aki_vals = lines[:,2].astype(float)
                lower_conf = np.vstack((lines[:,5],lines[:,6],lines[:,7])).T
                upper_conf  = np.vstack((lines[:,8],lines[:,9],lines[:,10])).T  


                if element_name not in atomh5.keys():
                    atomh5.create_group(element_name)
                    atomh5[element_name].attrs["Symbol"] = atom_dictionary[element_name]["symbol"]

                if ion_name not in atomh5[element_name].keys():
                    ion_group = atomh5[element_name].create_group(ion_name)

                    lambda0_set = ion_group.create_dataset('lambda0', data = lambda0)
                    lambda0_unit = self.SetUnit(vardescription = 'Rest frame wavelength [Å]',Lunit=1e-8,aFact=0.0,hFact=0.0)
                    for lambda_key in lambda0_unit.keys():
                        lambda0_set.attrs[lambda_key] = lambda0_unit[lambda_key]
                    lambda0_set.attrs['Source'] = 'NIST database'

                    fvalue_set  = ion_group.create_dataset('f-value', data = fvalue)
                    fvalue_unit = self.SetUnit(vardescription = 'Oscillator strength, dimensionless.',Lunit=1.0,aFact=0.0,hFact=0.0)
                    for fvalue_key in fvalue_unit.keys():
                        fvalue_set.attrs[fvalue_key] = fvalue_unit[fvalue_key]
                    fvalue_set.attrs['Source'] = 'NIST database'

                    aki_set = ion_group.create_dataset('Aki', data = aki_vals)
                    aki_unit = self.SetUnit(vardescription = 'Einstein Coefficients [s^-1]',Lunit=1.0,aFact=0.0,hFact=0.0)
                    for aki_key in aki_unit.keys():
                        aki_set.attrs[aki_key] = aki_unit[aki_key]
                    aki_set.attrs['Source'] = 'NIST database'

                    #dt = h5py.special_dtype(vlen=lower_conf.dtype)
                    lower_conf = ion_group.create_dataset('Lower Configuration', data = np.array(lower_conf,dtype='S'))
                    lower_conf.attrs['VarDescription'] = 'Lower level information. First column: Principal configuration. Second column: Principal term. Third column: J value'

                    upper_conf = ion_group.create_dataset('Upper Configuration', data = np.array(upper_conf,dtype='S'))
                    upper_conf.attrs['VarDescription'] = 'Upper level information. First column: Principal configuration. Second column: Principal term. Third column: J value'
                    
                    source_url = ion_group.create_dataset('url',data = line_url)
                    source_url.attrs['VarDescription'] = 'URL that was used to query the NIST database '
        atomh5.close()
        print("Done!")
    def add_line_to_hdf5(self,file_name="atom_dat.hdf5",element_name="Oxygen",ion_name="O VII",fvalue=0.696,lambda0=21.601690):
        """
        Adds a line for a specific ion to an HDF5 file.
        
        Args:
            file_name (str): Name of the HDF5 file.
            element_name (str): Name of the element.
            ion_name (str): Ion name in the format "ElementSymbol RomanNumeral".
            fvalue (float): Oscillator strength.
            lambda0 (float): Rest frame wavelength in Ångströms.
        """
        
        atomh5 = h5py.File(file_name,"a")
        if self.verbose:
            print("Adding "+ion_name+" with fval: "+str(fvalue)+" lambda0: "+str(lambda0)+" to: "+file_name)
        if element_name not in atomh5.keys():
            atomh5.create_group(element_name)
            
        if ion_name not in atomh5[element_name].keys():
            atomh5[element_name].create_group(ion_name)
            lambda0_set = atomh5[element_name][ion_name].create_dataset('lambda0', data = [lambda0])
            lambda0_unit = self.SetUnit(vardescription = 'Rest frame wavelength [Å]',Lunit=1e-8,aFact=0,hFact=0)

            for lambda_key in lambda0_unit.keys():
                lambda0_set.attrs[lambda_key] = lambda0_unit[lambda_key]
            lambda0_set.attrs['Source'] = 'Provided by the user'

            fvalue_set  = atomh5[element_name][ion_name].create_dataset('f-value', data = [fvalue])
            fvalue_unit = self.SetUnit(vardescription = 'Oscilator strength, dimensionless.',Lunit=1,aFact=0,hFact=0)

            for fvalue_key in fvalue_unit.keys():
                fvalue_set.attrs[fvalue_key] = fvalue_unit[fvalue_key]
            fvalue_set.attrs['Source'] = 'Provided by the user'        
            
        else:
            print("Oops, ion is already in the file")
                

        atomh5.close()


    def SetUnit(self,vardescription = 'text describing variable', Lunit=1, aFact=0, hFact=0):
        """
        Creates a dictionary containing metadata and conversion factors for a value.
        
        Args:
            vardescription (str): Description of the variable.
            Lunit (float): Conversion factor to CGS units.
            aFact (float): Exponent for the expansion factor (cosmological context).
            hFact (float): Exponent for the Hubble constant (cosmological context).
            
        Returns:
            dict: Dictionary with metadata and conversion factors.
        """

        
        return {'VarDescription': vardescription, 'CGSConversionFactor':Lunit, 'aexp-scale-exponent' :aFact, 'h-scale-exponent': hFact}



    def get_nist_line(self,ion2do = 'H I'):
        """
        Queries the NIST database for absorption line data of a specific ion.
        
        Args:
            ion2do (str): Ion name in the format "ElementSymbol RomanNumeral".
            
        Returns:
            tuple: (numpy.ndarray of lines, str of the query URL).
        """
        
        low_w = str(self.wavelength_low)
        upp_w = str(self.wavelength_high)
        grn_state = self.get_groundstate(ion2do)
        url = 'https://physics.nist.gov/cgi-bin/ASD/lines1.pl?spectra='+ion2do+'&limits_type=0&low_w='+low_w+'&upp_w='+upp_w+'&unit=0&de=0&I_scale_type=1&format=1&line_out=0&remove_js=on&no_spaces=on&en_unit=0&output=0&page_size=15&show_calc_wl=1&order_out=0&max_low_enrg=&show_av=3&max_upp_enrg=&tsb_value=0&min_str=&A_out=0&f_out=on&intens_out=on&max_str=&allowed_out=1&min_accur=&min_intens=&conf_out=on&term_out=on&J_out=on&submit=Retrieve+Data'
        respond = requests.get(url=url)

        self.check_for_query_errors(respond)

        table_post_pre = BeautifulSoup(respond.text, features='html5lib').find('pre')
        lines = table_post_pre.text.strip().split('\n')


        lines = self.format_lines_into_array(lines)
        lines = self.mask_ground_states(grn_state,lines)

        strong_line = None         
        intensities = np.array([float(self.clean_rel_ints(lines[:,1][i])) for i in range(len(lines))])

        if not all(intensities == 0.0):
            if self.verbose:
                print("Found most intense")
            most_intense_indx = np.argsort(intensities)[::-1][0]
            strong_line = lines[most_intense_indx]
        strong_line = None         

        wave_mask   = np.argsort(lines[:,0].astype('float'))[::-1]
        wavelengths = lines[:,0].astype('float')[wave_mask]


        ii=0
        lines_found = []
        delta_lambda =  self.delta_lambda 
        lines_wv_masked  = lines[wave_mask]

        filter_doubles = np.array([ s[-3].isdigit() for s in lines_wv_masked ])
        if ion2do == 'H I':

            lines_found = lines_wv_masked[filter_doubles]

        else:
            while ii<len(wavelengths):
                # we group the wavelengths that are within certain some ranges in this case 0.1
                mask_indx = np.where(abs((wavelengths[ii]-wavelengths))<delta_lambda)[0]
                nsublines = len(mask_indx)
                lines_wv_masked  = lines[wave_mask][mask_indx]
                # In the Nist data base sometimes we have two close wavelengths and the sum of them we identify them so we take the summed ones instead of them
                condition_array = np.array([lines_wv_masked[i][-3:]==["X","0","0"] for i in range(nsublines)]) 
                pattern_to_have = np.array([False,True,True])
                is_in_list = np.all(pattern_to_have == condition_array, axis=1)
                
                
                # If we find the summed line we take it and ignore the other lines. 
                if any(is_in_list):
            #        print("summed line found")
            #        print(lines[wave_mask][mask1_indx][is_in_list])
                    lines_found.append(lines_wv_masked[is_in_list])
                    ii += nsublines
                else:
                    if nsublines == 1:
                        lines_found.append(lines_wv_masked)
                    else:
                        fvalue_first = float(lines_wv_masked[0][3])
                        while len(lines_wv_masked)>1:
                            fvalue_first += float(lines_wv_masked[1][3])
                            lines_wv_masked = np.delete(lines_wv_masked,1,0)

                        lines_wv_masked[0][3] = str(fvalue_first)
                        lines_found.append(lines_wv_masked)
                    ii += nsublines                                          
            lines_found = np.array(lines_found)       

        fvals_float = np.vstack(lines_found)[:,3].astype('float')
        wvvals_float = np.vstack(lines_found)[:,0].astype('float')
        maxfval_indx = np.argsort(fvals_float*wvvals_float)[::-1]

        
        if lines_found.size != 0:       
        
            lines_fval_sorted = np.vstack(lines_found)[maxfval_indx]
        
        else:
            lines_fval_sorted = np.array([])
            
        if strong_line is None:
            strong_line = lines_fval_sorted[0]
        
        strong_indx = np.where(np.all(strong_line == lines_fval_sorted, axis=1))[0]
        weak_lines = lines_fval_sorted
        
        if strong_indx.size:
            weak_lines = np.delete(weak_lines,strong_indx,0)
        lines = np.vstack((strong_line,weak_lines))
        return lines,url