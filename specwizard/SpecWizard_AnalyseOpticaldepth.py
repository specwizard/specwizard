import h5py
import numpy as np
import random
from scipy.integrate import cumtrapz, simps
from scipy.interpolate import interp1d
from scipy.signal import convolve
import Phys

constants = Phys.ReadPhys()

specfile = 'spectra_part_los_z3.027.hdf5'
specdir  = 'path/for/specwizard/output.hdf5'

class Analyse_Opticaldepth:
    ''' Methods to analyse optical depth as a function of velocity '''
    def __init__(self,
                 specdir="./" ,
                 specfile=" ",
                 element='Hydrogen',
                 ion='H I',
                 wizard=None):
        
        self.dirname = specdir
        self.fname = specfile
        self.element = element
        self.ion = ion
        if wizard != None:
            self.nsight = wizard['sightline']['nsight']
            self.z = wizard['Header']['Cosmo']['Redshift']

            # determine number of pixels for each sight line
            #self.variable = 'LOS_0/' + self.element + '/' + self.ion + '/' + 'optical depth-weighted'
            #self.npix = hfile[self.variable].attrs.get("npix")
            self.header = wizard['Header']
            # determine pixel size
            self.pix_kms = wizard['pixkms']

            # determine length of spectrum
            self.box_kms = wizard['sightline']['Boxkms']['Value']
            self.npix    = int(self.box_kms/self.pix_kms) + 1 
            #self.box_kms = self.box_kms['Value']
            
                    # create velocity array
            self.vHubble = np.arange(self.npix) * self.pix_kms
            
            
        
        else:
            self.header = self.ReadHeader()

            print("Using output file")
            fname = self.dirname + '/' + self.fname
            hfile = h5py.File(fname, 'r')

            # determine number of sightlines
            self.nsight = hfile['Header'].attrs.get('number_of_sightlines')

            # determine redshift
            self.z = hfile['Header'].attrs.get('Redshift')

            # determine number of pixels for each sight line
            self.variable = 'LOS_0/' + self.element + '/' + self.ion + '/' + 'optical depth-weighted'
            self.npix = hfile[self.variable].attrs.get("npix")

            # determine pixel size
            self.pix_kms = hfile[self.variable + '/pixel_kms'][...]

            # determine length of spectrum
            self.box_kms = self.ReadVariable('LOS_0/Box_kms')
            self.box_kms = self.box_kms['Value']
                    # create velocity array
            self.vHubble = np.arange(self.npix) * self.pix_kms

    def ReadHeader(self):
        groupname = "Header"
        hfile = h5py.File(self.dirname + self.fname, "r")
        group = hfile[groupname]
        grp_dict = {}
        for k in group.attrs.keys():
            grp_dict[k] = group.attrs[k]
        hfile.close()
        header = dict(grp_dict)

        # read extra variables
        box = self.ReadVariable(varname="Header/Box")
        header["Box"] = box

        return header

    def ReadGroup(self, groupname='Header'):
        # read all entries for this particular hdf5 group
        hfile = h5py.File(self.dirname + self.fname, "r")
        group = hfile[groupname]
        grp_dict = {}
        for k in group.attrs.keys():
            grp_dict[k] = group.attrs[k]
        hfile.close()
        return dict(grp_dict)

    def ReadVariable(self, varname='Header/Box'):
        info = self.ReadGroup(groupname=varname)
        hfile = h5py.File(self.dirname + self.fname, "r")
        values = np.array(hfile[varname][...])
        hfile.close()
        return {'Value': values, 'Info': info}

    def read_optical_depths(self):
        """
        Read optical depths from an HDF5 file.

        Returns:
            dict: A dictionary containing information about the optical depths.

        """
        # Construct the full file path
        fname = self.dirname + '/' + self.fname

        # Open the HDF5 file
        hfile = h5py.File(fname, 'r')

        # Create an array to store results
        result = []
        var_path = f'LOS_{0}/{self.element}/{self.ion}/optical depth-weighted/'
        self.line_f0 = hfile[var_path+'f-value'][...]
        self.line_l0 = hfile[var_path+'lambda0'][...]
        # Iterate over sightlines and read optical depths
        for sight in np.arange(self.nsight):
            try:
                varname = f'LOS_{sight}/{self.element}/{self.ion}/optical depth-weighted/Optical depths'
                values = np.array(hfile[varname][...])
                result.append(values)
            except KeyError:
                print(f"Warning! Problem reading sightline #{sight}")

        # Close the HDF5 file
        hfile.close()

        # Create a dictionary with information about the optical depths
        result_dict = {
            'Info': "List of Optical depths (tau)",
            'Value': result,
        }

        return result_dict


    def mean_transmission(self, ODs, scale_factor=1.0):
        '''
        Calculate the mean transmission from optical depths.

        Args:
            ODs (numpy.ndarray or list): Optical depths. If a list, it should contain 1D arrays.
            scale_factor (float, optional): Scaling factor for optical depths. Defaults to 1.0.

        Returns:
            float: Mean transmission value.

        '''
        # Convert ODs to a numpy array
        ODs = np.array(ODs)

        # Check the dimensionality of ODs
        if ODs.ndim == 1:
            # If ODs is a 1D array, calculate the mean transmission
            return np.mean(np.exp(-scale_factor * ODs))
        else:
            # If ODs is a list of arrays, concatenate them and calculate the mean transmission
            AllODs = np.concatenate([ODs[i][:] for i in range(len(ODs))])
            return np.mean(np.exp(-scale_factor * AllODs))


    def scale_tau(self, ODs, meanflux=0.9, accuracy=1e-3, maxsteps=1500):
        '''
        Determine the scale factor for the optical depth so that mean transmission is equal to meanflux
        within the specified accuracy.

        Args:
            ODs (numpy.ndarray or list): Optical depths. If a list, it should contain 1D arrays.
            meanflux (float, optional): Target mean flux. Defaults to 0.9.
            accuracy (float, optional): Desired accuracy for mean transmission. Defaults to 1e-3.
            maxsteps (int, optional): Maximum number of steps for convergence. Defaults to 1500.

        Returns:
            float: Scale factor for the optical depth.

        '''
        # Initialize scale factor and calculate initial mean transmission
        scale_factor = 1.0
        meansim = self.mean_transmission(ODs, scale_factor=scale_factor)

        # Iteratively adjust the scale factor until desired accuracy is reached or maximum steps are reached
        nsteps = 0
        while (np.abs(meansim - meanflux) > accuracy) and (nsteps < maxsteps):
            nsteps += 1
            scale_factor *= np.sqrt(meansim / meanflux)
            meansim = self.mean_transmission(ODs, scale_factor=scale_factor)

        # Print a warning message if maximum steps are reached without reaching the desired accuracy
        if nsteps == maxsteps:
            print("Warning! Accuracy was not reached!")

        return scale_factor


    def BootStrapped_Flux(self,ODs,nsteps=5000):
        '''SimpleBootstrap procedure to get mean flux and standard deviation
           Input: 
           -Ods: List of optical depths
           -nsteps: number of iterations
           Output: Dictionary with:
           -Mean: Mean bootstraped flux 
           -STD : Standard Deviation
        '''
        bootstraped_flux = {}
        all_flux  = np.exp(-np.array(ODs))
        bootstraped_mean_flux = np.array([np.mean(flx) for flx in all_flux])
        hist_means = [np.mean([ random.choice(bootstraped_mean_flux) for i in range(len(bootstraped_mean_flux))]) for i in range(nsteps)]

        bootstraped_flux['Info']  = "Bootstrapped mean and standard deviation" 
        bootstraped_flux['Mean']  = np.mean(hist_means)#{}
        bootstraped_flux['STD']   = np.std(hist_means)#{}
        
        return bootstraped_flux
    def Onorbe_EffectiveTau(self, z):
        # fit from https://ui.adsabs.harvard.edu/abs/2017ApJ...837..106O/abstract
        A = 0.00126
        B = 3.294
        return A * np.exp(B * np.sqrt(z))
    
    def Kim2020_PowerLawTauFit(self,z,ver=2):
        ''' Power Laws presented in Kim2020, in the paper they present two power law fits. A simple one (ver=1) and a picewise (ver=2) 
            Input: z: Redshift, ver: (default:2) The type of the powerlaw from the paper.
            Output: tau
        '''
        if ver==1:
            A0    = -0.0060
            alpha = 2.87
        
        if ver==2:
            if z < 1.5:
                A0    = -0.0145
                alpha = 1.86
            if z >1.5:
                A0    = -0.0040
                alpha = 3.18
            
            
        return -1*A0*(1+z)**alpha
            
    def kim2007_TauFit(self,z):
        ''' Power Law presented in Kim2007 '''
        if (z>1.7) and (z<4):
            A0    = 0.0023
            alpha = 3.65 
            return A0*(1+z)**alpha
        

    def specwizard_tau_fit(self, z):
        """
        Fit the effective optical depth (tau)

        Args:
            z (float): Redshift value.

        Returns:
            float: Computed effective optical depth (tau).
        """
        # Coefficients for the fit function
        a0 = 0.04300355
        b0 = -0.03759746
        c0 = 0.01261603

        # Compute the effective optical depth (tau) using the fit function
        tau = a0 * (1 + z) + b0 * (1 + z)**2 + c0 * (1 + z)**3

        return tau

    def flux_ps(self, flux, meanf, N=None, V=None):
        """
        Calculate dimensionless flux power spectra.

        Args:
            flux (numpy.ndarray): 1D array representing the input flux data.
            meanf (float): Mean value of the flux data.
            N (int, optional): Number of pixels. Defaults to None, in which case it uses self.npix.
            V (float, optional): Box size in km/s. Defaults to None, in which case it uses self.box_kms.

        Returns:
            Tuple (numpy.ndarray, numpy.ndarray): Tuple containing the frequency values and corresponding dimensionless flux power spectra.

        """
        # Set default values if not provided
        if N is None:
            N = self.npix
        if V is None:
            V = self.box_kms

        # Calculate velocity bin size
        dv = V / N

        # Calculate frequency values
        freqs = np.fft.fftfreq(N)
        freqs *= (2 * np.pi / dv)

        # Sort frequencies and remove negative values
        indx = np.argsort(freqs)
        indx = indx[freqs[indx] >= 0]
        freqs = freqs[indx]

        # Calculate delta and perform Fourier transform
        delta = (flux - meanf) / meanf
        fourier = np.fft.fft(delta)[indx]
        fourier /= N

        # Calculate power spectrum
        pwr_spec = np.abs(fourier)**2
        pwr_spec *= V

        # Calculate dimensionless flux power spectra
        kPS = pwr_spec * freqs / np.pi

        return freqs, kPS

    
    def bin_PS(self,kphys,kPk):
        kbins     = np.arange(np.log10(0.0005), np.log10(0.07), 0.2)
        ksmall    = [0.001, 0.005, 0.007, 0.009] 
        klarge    = 10**np.arange(-2, -1, 0.1)
        kbins     = np.concatenate([ksmall, klarge])
        kcent     = 0.5*(kbins[1:]+kbins[0:-1])
        kPkbinned = np.zeros_like(kcent)
        for i in np.arange(len(kcent)):
            k1       = kbins[i]
            k2       = kbins[i+1]

            combined_kPk = []

            kphys    = kphys 
            kPk      = kPk
            mask     = ((kphys>=k1) & (kphys<k2))
            allkPk   = np.concatenate([kPk[i][mask] for i in range(len(kPk))])
            combined_kPk.append(allkPk)
            allkPk = np.concatenate([combined_kPk[i][:] for i in range(len(combined_kPk))])

            if allkPk.size > 0:
                kPkbinned[i] = np.mean(allkPk)
            else:
                kPkbinned[i] = 0 
        return kcent,kPkbinned    


    def even(self,n):
        # return true if input n is even, else returns false
        if (np.int(n/2)*2 == n):
            return True
        else:
            return False
    def Convolve(self,flux, FWHM):
        ''' Convolve the spectrum with instrumental broadining
         Input: 
           -spectrum: dictionary, containing
              -L:       linear extent of the array wavelength range
              -wave:    wavelength or Hubble velocity
              -flux:    flux
           -FWHM:       full-width at half maximum of the Gaussian line-spread function
        Output: the convolved spectrum in the form of a dictionary
        '''
        
#        wave = spectrum["wave"]
        wave = self.vHubble
        L    = wave[-1]
        dx   = wave[1:] - wave[0:-1]
#        dx   = np.concatenate([dx, [spectrum['L']-wave[-1]]])
        dx   = np.concatenate([dx, [L-wave[-1]]])

        wave = wave + 0.5 * dx
#        flux = spectrum["flux"]

        # create Gaussian
        sigmaG    = FWHM / (2*np.log(2))  # Gaussian's standard deviation
        n         = len(wave)
        if self.even(n):
            nc = np.int(n/2-1)
        else:
            nc = np.int((n+1)/2)-1
        waveG     = wave - wave[nc]
        Gauss     = np.exp(-(waveG)**2 / (2*sigmaG**2)) / (np.sqrt(2.*np.pi*sigmaG**2))
        #
        convolved = convolve(flux, Gauss, mode='same') / np.sum(Gauss)
        # convolved = flux
        #
        newspectrum = {}
        newspectrum["wave"] = wave
        newspectrum["flux"] = convolved
        newspectrum["FWHM"] = FWHM
        return newspectrum


    def Rebin(self,spectrum, wave):
        ''' rebin the function y(x) to the new bins xnew
         Interpolation is performed such that the mean of the function is conserved
         Input: 
           -spectrum: dictionary, containing
              -L:       linear extent of the wavelength range
              -wave:    wavelength or Hubble velocity
              -flux:    flux
           -wave:       new wavelenght or Hubble velocioty to rebin to 
           '''
        # determine pixel size
        L    = self.box_kms #spectrum["L"]
        npix = len(spectrum) #spectrum["npix"]
        pix  = L/npix
        x0   = np.arange(npix) * pix
        vHubble = self.vHubble
        dx = vHubble[1:] -vHubble[0:-1]
        dx = np.concatenate([dx, [L - x0[-1]]])
        #
        xR = np.copy(vHubble) + dx  # right side of pixel

        # cumulative sum of the function y0(x0)
        f = np.cumsum(spectrum * dx)

        # new pixels
        dxnew = wave[1:] - wave[0:-1]  # size of pixels
        dxnew = np.concatenate([dxnew, [L - wave[-1]]])
        xnewR = np.copy(wave) + dxnew  # right side of pixels

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

    def Noise(self,spectrum, SNR=100, seed=-1, Poisson=True):
        ''' generate Gaussian noise 
        Input:
           spectrum: spectrum dictionary
           SNR: value of the signal-to-noise at the continuum level
           seed: random seed
                seed<=0: randomly initialize random number generator
                seed>0: set seed
          Poisson: if True, generate Poisson-distributed noise with zero mean and this SNR
                   if False, generate Gaussian-distributed noise with zero mean and this SNR
        Output: for each entry in flux, return the corresponding value of the noise
        '''
        if seed <=0:
            np.random.seed(None)
        else:
            np.random.seed(seed)

        # compute standard deviation
        sigma = np.sqrt(1./SNR)

        # generate noise
        noise = None
        
    def add_instrument_gaussian(self,tau,FWHM_kms = 2.6,vx =0):
        '''
        Convolves (using scipy.signal.convolve) to the optical depth a gassian, simulating instrument noise.
        input: Optical depth, Full width maximum (default value 2.6 km/s)
        output: Convolved signal
        '''
        if len(vx) == 0:
            vx = self.vkms
            
        sigma    = FWHM_kms/(2*np.sqrt(2*np.log(2)))
        mu       = np.mean(vx)
        gaussian = np.exp(-np.power(vx - mu, 2.) / (2 * np.power(sigma, 2.)))
        tau_conv = convolve(tau,gaussian,mode='same') #/ np.sum(gaussian)
        
            

        return tau_conv
    def read_out_noise(self,rebinned_flux, variance=0.05,seed=-1):
        '''
        Takes the rebinned flux and add the noise coused by the read-out noise with certain variance.
        input: Rebinned flux, variance 
        '''
        if seed != -1:
            np.random.seed(seed)
            
        sigma = np.sqrt(variance)
        
        noise = np.random.normal(loc=0,scale = sigma, size = len(rebinned_flux))
        
        return noise
    
    def sh_noise(self,rebinned_flux, mean=1,seed=-1):
        '''
        Takes the rebinned flux and add the noise coused by the read-out noise with certain variance.
        input: Rebinned flux, variance 
        '''
        if seed != -1:
            np.random.seed(seed)
      
        noise = np.random.poisson(lam=mean, size = len(rebinned_flux))
        
        return noise
        

    def gauss_decomposition(self,  ODs,  tresh=[14, 15, 16]):
        """Perform Gaussian decomposition on spectral lines.

        Args:
            ODs: List of optical depths.
            tresh: List of thresholds for selecting data.

        Returns:
            Tuple: Dictionaries containing optical depth and flux results.
        """

        # Initialize threshold values for column density
        cds_tresh = 10**np.array(tresh)

        # Dictionaries to store optical depth and flux results for each threshold
        tau_f = {}
        flux_f = {}
        sigma = self.get_sigma()
        # Initialize entries in the dictionaries for each threshold
        for trsh in tresh:
            tau_f[trsh] = []
            flux_f[trsh] = []
        tau_f['total'] = []
        flux_f['total'] = []
        pixkms = self.pix_kms
        boxkms = self.box_kms

        # Loop over each spectrum in the input list ODs
        ODs = ODs['Value']
        for i in range(len(ODs)):
            # Extract necessary information about the spectrum
            optical_depth = ODs[i]
            npix = len(optical_depth)
            velocity = np.arange(npix) * pixkms

            # Decompose the spectral line
            line_maxima, line_start, line_end, column, column2, equv_w, nHtot = self.line_decomp(
                velocity, optical_depth, boxkms, sigma)

            # Create an array to store column density values
            cdens = np.zeros_like(velocity)

            # Populate the cdens array with column density values from the decomposition
            for start, end, cd in zip(line_start, line_end, column2):
                cdens[start:end] = cd

            # Loop over each column density threshold
            for cds in cds_tresh:
                key_trsh = np.log10(cds)

                # Create a mask for pixels below the current column density threshold
                mask = cdens < cds

                # If any pixels satisfy the condition, calculate optical depth and flux
                if any(mask):
                    tau_f[key_trsh].append(-np.log(np.mean(np.exp(-optical_depth[mask]))))
                    flux_f[key_trsh].append(np.exp(-optical_depth[mask]))
                else:
                    # If no pixels satisfy the condition, set the optical depth to 0
                    tau_f[key_trsh].append(0)

            # Calculate and append total optical depth and flux for the entire spectrum
            tau_f['total'].append(-np.log(np.mean(np.exp(-optical_depth))))
            flux_f['total'].append(np.exp(-optical_depth))

        # Finalize the results
        for cds in cds_tresh:
            key_trsh = np.log10(cds)
            # Recalculate the total optical depth for each threshold
            tau_f[key_trsh] = -np.log(np.mean(np.exp(-np.array(tau_f[key_trsh]))))

        # Recalculate the total optical depth for all spectra
        tau_f['total'] = -np.log(np.mean(np.exp(-np.array(tau_f['total']))))

        # Return the dictionaries containing the final results
        return tau_f, flux_f

    def line_decomp(self, velocity, optical_depth, boxkms, sigma):
        """Perform line decomposition for a spectral line.

        Args:
            velocity (array): Array of velocities.
            optical_depth (array): Array of optical depths.
            boxkms (float): Width of the spectral box in kilometers per second.
            sigma (float): Parameter used in calculations.

        Returns:
            tuple: Tuple containing various quantities related to the identified spectral lines.

        """
        # Calculate the number of pixels
        npix = len(velocity)
        
        # Calculate the velocity resolution per pixel
        pixkms = velocity[1] - velocity[0]
        
        # Calculate the total hydrogen column density
        nHtot = np.sum(optical_depth) * (pixkms * 1e5) / sigma

        # Extend arrays for analysis
        vel = np.concatenate((velocity - boxkms, velocity, velocity + boxkms))
        tau = np.concatenate((optical_depth, optical_depth, optical_depth))

        # Find minima and maxima in the optical depth
        minima = np.where(((tau <= np.roll(tau, -1)) & (tau < np.roll(tau, +1))) |
                        ((tau < np.roll(tau, -1)) & (tau <= np.roll(tau, +1))))[0]
        maxima = np.where((tau > np.roll(tau, -1)) & (tau > np.roll(tau, +1)))[0]

        # Handle cases with no maxima
        if len(maxima) == 0:
            return [[], [], [], [], [], [], 0]

        # Adjust indices if the first maximum occurs before the first minimum
        if maxima[0] < minima[0]:
            maxima = maxima[1:]
            minima = minima[:-1]

        # Define line start and end indices
        line_start = minima[0:-1]
        line_end = minima[1:]

        # Calculate derived quantities
        dderiv = (np.roll(tau, +1) + np.roll(tau, -1) - 2 * tau) / (pixkms * 1e5) ** 2
        vel_zero = vel[maxima]
        tau_zero = tau[maxima]
        bpar2 = -2 * tau_zero / dderiv[maxima]
        column = tau_zero * np.sqrt(np.pi * bpar2) / sigma

        # Calculate column density using a different method
        ctau = np.cumsum(tau) * (pixkms * 1e5)
        column2 = ctau[minima[1:]] - ctau[minima[:-1]]
        column2 /= sigma

        # Calculate equivalent width
        w = 1 - np.exp(-tau)
        cw = np.cumsum(w) * (pixkms)
        equv_w = cw[minima[1:]] - cw[minima[:-1]]

        # Mask velocities within a specified range
        vel_mask = np.where((vel_zero >= 0) & (vel_zero < boxkms))[0]

        # Extract relevant quantities based on the velocity mask
        line_maxima = maxima[vel_mask] - npix
        line_start = line_start[vel_mask] - npix
        line_end = line_end[vel_mask] - npix
        column = column[vel_mask]
        column2 = column2[vel_mask]
        equv_w = equv_w[vel_mask]

        # Return the results
        return line_maxima, line_start, line_end, column, column2, equv_w, nHtot

    def get_sigma(self):
        """gets sigma

        Args:
            None
        Returns:
            tuple: Tuple containing various quantities related to the identified spectral lines.

        """
        line_f     = self.line_f0
        line_l0    = self.line_l0
        sigmaT     = constants['sigmaT']

        sigma_line   = np.sqrt(3*np.pi*sigmaT/8) * line_f * line_l0 * 1e-8 * constants['c']
        return sigma_line