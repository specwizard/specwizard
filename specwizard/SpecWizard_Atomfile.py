import numpy as np
import requests
from bs4 import BeautifulSoup
import h5py
import traceback
import mendeleev

class Atomfile:
    """ This Class creates the atom file that contains the atomic transition data that will be use to create the spectra line of ions.
        The data is retrieved from the NIST line database. https://physics.nist.gov/PhysRefData/ASD/lines_form.html 

    """


    def __init__(self):

        self.nist_url = 'https://physics.nist.gov/cgi-bin/ASD/lines1.pl?'

        self.nist_grnd_url = 'https://physics.nist.gov/cgi-bin/ASD/ie.pl?'

    # we generate an atom direcotry with information about the element and ions

        romnumerals    = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 
                    'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI', 'XVII',
                'XVIII', 'XIX', 'XX', 'XXI', 'XXII', 'XXIII', 'XXIV', 'XXV', 'XXVI']
        
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
                atom_dictionary[element]["lines"].append(element_short[i]+" "+romnumerals[j])

        self.atom_dictionary = atom_dictionary 

    def format_line(self,line):
        """
        This function add's a zero whenever a column from the table has an empty space. 
        """
        for i,ln in enumerate(line):
            if len(ln.strip()) == 0:
                line[i] = '0'
            else:
                line[i] = ln.strip()
        return line    
    
    def clean_rel_ints(self,value):
        """
        In the NIST database the value of intensity is often found accompanied with certain characters. This function clean those characters so it can be converted into a string
        """
        rel_ints_markers = ['*',':','-',',','a','bl','b','B','c','d','D','E','N','f','F','e','g','G','H','hfs','h','i','j','l','m','o','O','p','q','r','s','S','I','t','u','w','x','(',')','IV/2',' ','V/2','?']
        for marker in rel_ints_markers:
            value = value.replace(marker,"")
        
        if value == '':
            value = '0'
        return value    
    
    def format_lines_into_array(self,lines):
        '''
        This function takes the string lines from the query, and formats it into a numpy array.
        '''
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
        Checks for possible errors from the query response. It can be that the line is not available in the current wavelength range or that the user give a wrong ion key e.g (HI,h I,He IX)
        """
        if 'No lines are available in ASD with the parameters selected' in response.text:
            
            raise Exception("Line not in the NIST database with current settings.")
            
        if 'Unrecognized token.' in response.text:
            
            raise Exception("Unrecognized ion")

    def mask_ground_states(self, ground_state,lines):
        """
        Will check for the lines in the ground state, and only keep those. In case of not finding any, it will raise an exception. 
        """
    
        ground_mask =lines[:,5]==ground_state

        if any(ground_mask):
            return lines[ground_mask]
        else:
            raise Exception("Ground state not found!")



    def get_groundstate(self,ion_to_q):
        """
        This Function query the nist data base to get the ground level configuration, this in order to only get the lines in such configuration.

        Args:
            ion_to_q (str): Ion to query, it must be in roman numerals and with space e,g "H I", "Si IV" 

        Returns:
            str : String of the ground level configuration. 
        """

        url_grnd = self.nist_grnd_url+'spectra='+ion_to_q+'&units=1&format=1&remove_js=on&order=0&conf_out=on&e_out=0&submit=Retrieve+Data' 
        respond = requests.get(url=url_grnd)
        table_post_pre = BeautifulSoup(respond.text, features='html5lib').find('pre')
        lines = table_post_pre.text.strip().split('\n')
        ground_state = lines[3].replace('|','').split()[0]
        
        return ground_state        

    def create_hdf5_from_nist(self,file_name="atom_dat.hdf5", wavelength_low_lim=200.0,wavelength_upper_lim=8000.0,wavelength_resolution=0.1):
        """
        This functions writes a hdf5 file containing the atomic transition of the strongest lines of each Ion. The information is query from the NIST data base. 
        The default elements it queries are Hydrogen, Deuterium, Helium, Carbon, Oxygen, Silicon, Iron, Magnesium

        Args:
            file_name (str): file_name.hdf5 
            wavelength_low_lim (float): Lower wavelength limit from the spectrograph, this will delimit the lines that we query in the database. In Amstrongs [Å]
            wavelength_upper_lim (float): Upper wavelength limit from the spectrograph, this will delimit the lines that we query in the database. In Amstrongs [Å]
            wavelength_resolution (float):  Resolution of the spectrograph. This translates into, if there are lines that are within this range, we will take it as one summing their f-value. In Amstrongs [Å]
        """
        atom_dictionary = self.atom_dictionary

        self.wavelength_low = wavelength_low_lim
        self.wavelength_high = wavelength_upper_lim
        self.delta_lambda    = wavelength_resolution
        atomh5 = h5py.File(file_name,"a")

        for element_name in atom_dictionary.keys():
            for ion_name in atom_dictionary[element_name]["lines"]:

                print("Writting "+element_name+ " "+ ion_name+ "...")
                try:
                    lines, line_url= self.get_nist_line(ion2do=ion_name)
                except Exception as exc:

                    #print(traceback.format_exc())
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
    
    def add_line_to_hdf5(self,file_name="atom_dat.hdf5",element_name="Oxygen",ion_name="O VII",fvalue=0.696,lambda0=21.601690):
        """
        This function adds to (or create) to a hdf5 file a line for a particular ion

        Args:
            file_name (str): file name to add line
            element_name (str): name of element
            ion_name (str): Ion name in the format of abreviation+*space*+roman_numeral e.g H I, Fe VII
            fvalue (float): Osilator strenght of the ion. 
            lambda0 (float): Rest frame wavelength in Amstrongs (Å)
        """
        
        atomh5 = h5py.File(file_name,"a")
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
        '''
        Creates dictionary to add context and conversion factors to a value.

        Args:
            vardescription (str): Text description of the variable
            Lunit (float): Conversion value for CGS
            aFact (float): expansion factor exponent. Needed for some cosmological simulations. 
            hFact (float): hubble constant exponent. Needed for some cosmological simulations. 

        Returns:
            dict : Python dictionary with the formated values. 

        '''
        return {'VarDescription': vardescription, 'CGSConversionFactor':Lunit, 'aexp-scale-exponent' :aFact, 'h-scale-exponent': hFact}



    def get_nist_line(self,ion2do = 'H I'):
        '''
        This function query the NIST database to obtain the ion absorption line information.

        Args:
            ion2do (str): Ion name as formated with the element chemical symbol+*space*+roman numeral e.g (H I, Fe XI)

        '''
        
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
            print("Found most intense")
            most_intense_indx = np.argsort(intensities)[::-1][0]
            strong_line = lines[most_intense_indx]

        wave_mask   = np.argsort(lines[:,0].astype('float'))[::-1]
        wavelengths = lines[:,0].astype('float')[wave_mask]



        ii=0
        lines_found = []
        delta_lambda =  self.delta_lambda 
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

        fvals_float = np.vstack(lines_found)[:,3].astype('float')                
        maxfval_indx = np.argsort(fvals_float)[::-1]
        if lines_found:
        
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