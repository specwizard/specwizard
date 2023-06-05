import h5py
import numpy as np
import random
from scipy.integrate import cumtrapz, simps
from scipy.interpolate import interp1d
from scipy.signal import convolve

specfile = 'spectra_part_los_z3.027.hdf5'
specdir  = 'path/for/specwizard/output.hdf5'

class Analyse_Opticaldepth:
    ''' Methods to analyse optical depth as a function of velocity '''
    def __init__(self,
                 specdir=" " ,
                 specfile=" ",
                 element='Hydrogen',
                 ion='H I',
                 Wizard=None):
        
        self.dirname = specdir
        self.fname = specfile
        self.element = element
        self.ion = ion
        
        if Wizard != None:
            self.nsight = Wizard['sightline']['nsight']
            self.z = Wizard['Header']['Cosmo']['Redshift']

            # determine number of pixels for each sight line
            #self.variable = 'LOS_0/' + self.element + '/' + self.ion + '/' + 'optical depth-weighted'
            #self.npix = hfile[self.variable].attrs.get("npix")
            self.header = Wizard['Header']
            # determine pixel size
            self.pix_kms = Wizard['pixkms']

            # determine length of spectrum
            self.box_kms = Wizard['sightline']['Boxkms']['Value']
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

    def Read_Opticaldepth(self):

        fname = self.dirname + '/' + self.fname
        hfile = h5py.File(fname, 'r')

        # create array to store results
        result = []
        for sight in np.arange(self.nsight):
            varname = 'LOS_{}'.format(
                sight
            ) + '/' + self.element + '/' + self.ion + '/' + 'optical depth-weighted/Optical depths'
            values = np.array(hfile[varname][...])
            result.append(values)

        # create velocity array

        hfile.close()

        result = {
            'Info': "List of Optical depths (tau)",
            'Value': result,
                 }
        return result

    def Mean_Transmission(self, ODs, scale_factor=1.0):
        ''' determine mean optical depth '''

        AllODs = np.concatenate([ODs[i][:] for i in range(len(ODs))])
        return np.mean(np.exp(-scale_factor * AllODs))

    def Scale_Tau(self, ODs, meanflux=0.9, accuracy=1e-3,maxsteps = 1500):
        # Determine the scaline factor for the optical depth so
        # that mean transmission is equal to meanflux within the specified accuracy
        scale_factor = 1.0
        meansim = self.Mean_Transmission(ODs, scale_factor=scale_factor)
        nsteps = 0
        while (np.abs(meansim - meanflux) > accuracy) and (nsteps < maxsteps):
            nsteps += 1
            scale_factor *= np.sqrt(meansim / meanflux)
            meansim = self.Mean_Transmission(ODs, scale_factor=scale_factor)
            
        if nsteps ==maxsteps:
            print("Warning! accuracy was not reached!") 
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
    
    def FluxPS(self,flux,meanf, N=None, V=None):
        ''' Function that calculates the flux power spectrum
            Input: Flux (exp(-tau))
                   MeanFlux: average flux of all your spectra NOT of the individual spectra, or theoretical Flux eg from oñorbe fit
                   N: Number of pixels
                   V: velocity of the box
            output: Kvals in
                  : K*PS

        '''
        if N==None:
            N = self.npix
        if V==None:
            V = self.box_kms
        dv       = V / N
        freqs    = np.fft.fftfreq(N) 
        freqs   *= (2 * np.pi / dv)
        indx     = np.argsort(freqs)
        indx     = indx[freqs[indx] >= 0]
        freqs    = freqs[indx]

        delta    = (flux - meanf) / meanf 
        fourier  = np.fft.fft(delta)[indx]
        fourier  /= N
        pwr_spec = np.abs(fourier)**2
        pwr_spec *= V


        kPS      = pwr_spec * freqs 

        return freqs,kPS
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
    def Convolve(self,spectrum, FWHM):
        ''' Convolve the spectrum with instrumental broadining
         Input: 
           -spectrum: dictionary, containing
              -L:       linear extent of the array wavelength range
              -wave:    wavelength or Hubble velocity
              -flux:    flux
           -FWHM:       full-width at half maximum of the Gaussian line-spread function
        Output: the convolved spectrum in the form of a dictionary
        '''
        wave = spectrum["wave"]
        dx   = wave[1:] - wave[0:-1]
        dx   = np.concatenate([dx, [spectrum['L']-wave[-1]]])
        wave = wave + 0.5 * dx
        flux = spectrum["flux"]

        # create Gaussian
        sigmaG    = FWHM / (2*np.log(2))  # Gaussian's standard deviation
        n         = len(wave)
        if even(n):
            nc = np.int(n/2-1)
        else:
            nc = np.int((n+1)/2)-1
        waveG     = wave - wave[nc]
        Gauss     = np.exp(-(waveG)**2 / (2*sigmaG**2)) / (np.sqrt(2.*np.pi*sigmaG**2))
        #
        convolved = convolve(flux, Gauss, mode='same') / np.sum(Gauss)
        # convolved = flux
        #
        newspectrum = spectrum.copy()
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
        