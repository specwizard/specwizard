# %load SpecWizard_Lines_tom.py
import numpy as np
import unyt
import scipy.interpolate as interpolate
from scipy.signal import convolve
from scipy.special import erf
from astropy.modeling.functional_models import Voigt1D
from scipy.special import voigt_profile as VoigtSciPy
from .Phys import ReadPhys
constants = ReadPhys()
from unyt.dimensions import length, time, mass, temperature
from unyt import accepts

# convolve(in1, in2, mode='full', method='auto')[source]
class Lines:
    ''' Methods to compute optical depth as a function of velocity for a single absorber,
        implementing:
        - a Gaussian profile
        - a Voigt profile
    '''

    def __init__(self, v_kms =0.0, box_kms=-1.0, constants=constants, lambda0_AA=1215.67, f_value=0.4164,
                 naturalwidth_kms=6.06076e-3,verbose=False, periodic=True):

        self.constants    = constants
        self.verbose      = verbose
        self.v_kms        = v_kms     # velocity bins to evaluate optical depth [km/s]
        self.pix_kms      = v_kms[1] - v_kms[0] # pixel size
        self.periodic     = periodic
        self.box_kms      = box_kms             # velocity extent of spectrum
        self.npix         = len(self.v_kms)
        self.lambda0      = lambda0_AA * 1e-8  * unyt.cm # rest-wavelength           [cm]
        self.f_value      = f_value            # oscillator strength       [dimensionless]
        self.naturalwidth = naturalwidth_kms    # natural line width        [km/s]self.
        self.sigma        = self.constants["c"] * np.sqrt(3*np.pi*self.constants["sigmaT"] /8.) * self.f_value * self.lambda0
    def errfunc(self):
        # tabulated error function
        pix     = 1e-2
        vmax    = 100.
        npix    = int(vmax/pix)
        pix     = vmax / npix
        v       = np.arange(-npix, npix) * pix
        err     = 0.5 * (1 + erf(v))
        return v, err
    
    def IDinterpol(self, x, xp, fp, cumulative=True):
        ''' Interpolate the function fx(xp) to the points x, conserving the integral of fp.
            This assumes that both x and xp are equally spaced '''
        # extend x axis by one element
        dx   = x[-1] - x[-2]
        xnew = np.concatenate((x, [x[-1]+dx]))


        # compute culmulative sum of fp's and extend by one element
        if not cumulative:
            Fc   = np.concatenate(([0], np.cumsum(fp)))
        else:
            Fc   = np.concatenate(([0], fp))
        dX   = xp[-1]-xp[-2]
        Xc   = np.concatenate((xp, [xp[-1]+dX]))

        # interpolate cumulative sum
        fcnew = np.interp(xnew, Xc, Fc)
        # difference
        fnew  = (np.roll(fcnew, -1) - fcnew)[:-1]
        fnew  = fnew
        return fnew    

    @accepts(column_densities=length**-2,b_kms=length/time,vion_kms=length/time,Tions=temperature)
    def gaussian(self, column_densities = 0, b_kms = 0,vion_kms=0,Tions= 0, mass_densities=0):
        '''
            column_density = ion column density in ions / cm^2
            density        = mass density of the gas
            b_kms          = line-width parameter
            vion_kms       = peculiar velocity
            Tions          = temperature of the gas
            mass_densities = (total) density of the gas 
        '''

        naturalwidth_kms = self.naturalwidth    # natural line width        [km/s]
        f_value      = self.f_value
        lambda0      = self.lambda0  # rest-wavelength           [cm]
        pix_kms      = self.pix_kms       # pixel size
        periodic     = self.periodic

        # line cross section times speed of light, c [cm^2 * cm/s]
        sigma        = self.sigma.in_cgs().value
        
        # generate normalized error function
        verf, erf = self.errfunc()
        
        # extent velocity range
        pixel_velocity_kms = np.concatenate((self.v_kms - self.box_kms, self.v_kms, self.v_kms + self.box_kms))
        tau          = np.zeros_like(pixel_velocity_kms)
        densities    = np.zeros_like(pixel_velocity_kms) #* column_densities
        velocities   = np.zeros_like(pixel_velocity_kms) #* vion_kms
        temperatures = np.zeros_like(pixel_velocity_kms) #* Tions.units

        no_unyt_pixel_velocity_kms = pixel_velocity_kms.in_cgs().value

        #strip units for performance 
        col_dens_unit, no_unyt_column_densities = column_densities.in_cgs().units, column_densities.in_cgs().value
        b_unit, no_unyt_b_kms          = b_kms.in_cgs().units, b_kms.in_cgs().value
        vion_unit, no_unyt_vion_kms    = vion_kms.in_cgs().units, vion_kms.in_cgs().value
        T_unit, no_unyt_Tions          = Tions.in_cgs().units, Tions.in_cgs().value
        d_unit, no_unyt_mass_densities = mass_densities.in_cgs().units, mass_densities.in_cgs().value
        for column_density, b, vel, Tion, mass_density in zip(no_unyt_column_densities, no_unyt_b_kms, no_unyt_vion_kms,no_unyt_Tions, no_unyt_mass_densities):
            if column_density >0:
                # scale b-parameter
                v_line = b * verf

                # interpolate, and convert velocity from km/s to cm/s

                g_int   = column_density * sigma * self.IDinterpol(no_unyt_pixel_velocity_kms - vel, v_line, erf, cumulative=True) 
                # add
                tau          += g_int
                densities    += g_int * mass_density
                velocities   += g_int * vel
                temperatures += g_int * Tion            

        # normalize to pixel size
        pix_cms       = self.pix_kms.in_cgs().value
        tau          /= pix_cms
        velocities   /= pix_cms
        temperatures /= pix_cms
        densities    /= pix_cms
        
        # Give back units 
        densities    *= d_unit 
        velocities   *= vion_unit 
        temperatures *= T_unit
        nint = self.npix
        
        if periodic:
            tau          = tau[0:nint] + tau[nint:2*nint] + tau[2*nint:3*nint]
            pixel_velocity_kms = pixel_velocity_kms[nint:2*nint] 
            densities    = densities[0:nint] + densities[nint:2*nint] + densities[2*nint:3*nint]
            velocities   = velocities[0:nint] + velocities[nint:2*nint] + velocities[2*nint:3*nint]
            temperatures = temperatures[0:nint] + temperatures[nint:2*nint] + temperatures[2*nint:3*nint]
            
        else:
            tau   = tau[nint:2*nint]
            pixel_velocity_kms = pixel_velocity_kms[nint:2*nint] 
            densities    = densities[nint:2*nint] 
            velocities   = velocities[nint:2*nint]
            temperatures = temperatures[nint:2*nint]            
        mask = tau > 0
        
        #Normalize optical depth-weighted  quantities 
        densities[mask]     /=  tau[mask]
        velocities[mask]    /=  tau[mask]
        temperatures[mask]  /=  tau[mask]
        

        # compute total column density
        nh_tot = np.cumsum(tau)[-1] * pix_cms / self.sigma.in_cgs()

        # Set units for velocities, temperatures and densities
        vunit        = self.SetUnit(vardescription="Tau weighted ion peculiar velocity", aFact=0.0, hFact=0.0)
        velocities   = {'Value': velocities, 'Info': vunit}
        tunit        = self.SetUnit(vardescription="Tau weighted ion temperature", aFact=0.0, hFact=0.0)
        temperatures = {'Value': temperatures, 'Info': tunit}
        dunit        = self.SetUnit(vardescription="Tau weighted mass density", aFact=0.0, hFact=0.0)
        densities    = {'Value': densities, 'Info': dunit}
        tauunit      = self.SetUnit(vardescription="Ionic optical depth", aFact=0.0, hFact=0.0)
        tau          = {'Value': tau, 'Info': tauunit}

        
        spectrum = {'pixel_velocity_kms':pixel_velocity_kms,
            'optical_depth':tau,
            'optical_depth_densities':densities,
            'optical_depth_velocities':velocities,
            'optical_depth_temperatures':temperatures,
            'total_column_density':nh_tot}
        
        return spectrum
    

    def directgauss(self, column_density = 0, b_kms = 0, lambda0_AA=1215.67, f_value=0.4164, naturalwidth_kms=6.06076e-3,periodic=True):
        ''' Direct Gaussian line profile evaluated at centre of each pixel '''

        # cross section [cm^3/s]
        sigma     = self.constants["c"].value * np.sqrt(3*np.pi*self.constants["sigmaT"].value/8.) * f_value * lambda0 # cross section cm^3/s        
        
        # extent velocity range
        velocity_kms = np.concatenate((self.v_kms - self.box_kms, self.v_kms, self.v_kms + self.box_kms))
        tau          = np.zeros_like(velocity_kms)
        
        #
        pix2         = 0.5*self.pix_kms  # 1/2 of the size of a pixel
        for norm, b, vel in zip(column_density * sigma, b_kms, self.v_kms):
            
            # scale b-parameter
            if norm > 0:
                # By convention, velocity_kms is the start of the pixel. 
                # We evaluate the optical depth at the centre of the pixel
                dv = vel - (velocity_kms + pix2)   

                gauss     = 1./np.sqrt(np.pi*(1e5*b)**2) * np.exp(-(dv)**2/b**2)
                gauss    *= norm
                
                #
                tau      += gauss

        #
                
        # apply periodic boundary conditions if required
        nint = self.npix
        if periodic:
            tau          = tau[0:nint] + tau[nint:2*nint] + tau[2*nint:3*nint]
            velocity_kms = velocity_kms[nint:2*nint]
            
        # compute total column density
        nh_tot = np.cumsum(tau)[-1] * 1.e5 * self.pix_kms / sigma
        
        return velocity_kms, tau, nh_tot

    
    def convolvelorentz(self, phi):
            ''' return convolution of input line profile with Lorentzian
            input: phi(v_z): line profile function as a function of velocity, v_z
            output: phi(v_z) input line profile now convolved with natural line profile
            '''
            # The core of the Lorentzian needs to be sampled to 0.01 km/s to get an accurate
            # convolution. Therefore we interpolate the original profile to a finer velocity grid
            # before performing the convolution.

            # Create velocity bins for interpolation
            dv         = 1e-3         # pixel size in km/s
            vmin       = np.min(self.v_kms)
            vmax       = np.max(self.v_kms)
            nbins      = np.array((vmax-vmin)/dv,dtype=int)
            dv         = (vmax-vmin)/float(nbins)
            v_convolve = vmin + np.arange(nbins) * dv 

            

            # Create lorentz profile
            phi_fine    = np.interp(v_convolve, self.v_kms, phi)
            width       = self.naturalwidth
            v_bins      = v_convolve - np.mean(v_convolve)              # centre Lorenz at the centre of the velocity interval
            lorentz     = (1./np.pi) * width / (v_bins**2 + width**2)   
            lorentz     = lorentz                                 # convert to units of [s/cm]

            # The integral over the Lorentzian is normalized to unity, so that
            # sum(lorentz) dpix = 1, or the pixel size = 1/np.sum(lorentz)
            phi_fine    = convolve(phi_fine, lorentz, mode='same') / np.sum(lorentz)

            #
            result      = np.interp(self.v_kms, v_convolve, phi_fine)

            return result    #
        
    def SciPyVoigt(self, b_kms=10., v0_kms=0.0, lambda0_AA=1215.67, f_value=0.4164, naturalwidth_kms=6.06076e-3, periodic=True):
        ''' 
        return Voigt line-profile function, Equation 5
        this version uses the SciPy implementation of the Voigt function
        input:
             b_kms (float):     b-parameter [km/s]
             v0_kms (float):    velocity at line centre [km/s]
        output: 
             line profile (float or array) : line shape with unit column density [s/cm]
        
        '''
        # extend velocity array
        velocity_kms = np.concatenate((self.v_kms - self.box_kms, self.v_kms, self.v_kms + self.box_kms))
        u            = (velocity_kms-v0_kms) / b_kms
        
        # 
        sigma_G     = 1.0 / np.sqrt(2.0)           # variance of Gaussian
        gamma_L     = naturalwidth_kms / b_kms     # Half-width half-maximum of Lorenzian - the parameter "a"
        
        # evaluate Voigt profile
        vnorm       = VoigtSciPy(u, sigma_G, gamma_L, out=None)
        phi         = vnorm / b_kms                # SciPy returns a normalized Voigt profile, which includes the 1/sqrt(pi)

        # impose periodic boundary conditions
        nint = self.npix        
        if periodic:
            phi          = phi[0:nint] + phi[nint:2*nint] + phi[2*nint:3*nint]
            velocity_kms = velocity_kms[nint:2*nint]
        
        return velocity_kms, phi/1e5  # convert to units of [s/cm]
    
    
    
    
    def sigmaHI(self,hnu_eV=13.6):
        ''' Fit from Verner et al ('96) fit to the photo-ionization cross section
        Input: energy of the photon in eV
        Output: photo-ionization cross section in cm^2 '''

        barn   = 1e-24
        sigma0 = 5.475e4 * 1e6 * barn
        E0     = 0.4298
        xa     = 32.88
        P      = 2.963
        #
        energy = np.array(hnu_eV)
        x      = energy/E0
        sigma  = sigma0 * (x-1)**2 * x**(0.5*P-5.5) * (1.+np.sqrt(x/xa))**(-P)
        if isinstance(sigma, (list, tuple, np.ndarray)):
            mask        = energy < 13.6
            sigma[mask] = 0
        else:
            if energy < 13.6:
                sigma = 0.0    
        return sigma

    def convolveLymanLimit(self,tau_Lymanalpha):
        ''' Return Lyman limit optical depth corresponding to input Lyman-alpha optical depth '''
        vel_kms    = self.v_kms           # pixel velocities [km/s]
        pix_kms    = self.pix_kms    # pixel width [km/s]
        npix       = self.npix       # number of pixels
        constants  = self.constants
        lambda0    = 1215.67e-8      # Lya wavelength [cm]
        f_value    = 0.4164          # Lya f-value

        # compute Lya cross section [cm^2 km/s]
        sigma_a    = np.sqrt(3*np.pi*constants["sigmaT"].value/8.) * f_value * lambda0 * (constants["c"].value) / 1e5

        # generate velocity grid for convolution

        # use finer pixels
        pix     = 0.1 # pixel size in km/s
        npix    = int(np.max(vel_kms) / pix)
        pix     = np.max(vel_kms) / npix
        vel     = np.arange(-npix, npix) * pix
        tau     = np.interp(vel, vel_kms, tau_Lymanalpha)
        #hnu_eV  = 13.6 * (1. - vel * 1e5 / constants["c"])
        hnu_eV  = 13.6 * np.exp(-((vel*1e5)/constants["c"].value))
        sigma   = self.sigmaHI(hnu_eV=hnu_eV)
        tau_LL  = convolve(tau/sigma_a, sigma, mode='same')
        tau_LL *= pix

        #
        result  = np.interp(vel_kms, vel, tau_LL)

        #
        return result
        
    def SetUnit(self, vardescription='text describing variable',  aFact=1.0, hFact=1.0):
        ''' 
        Set the unit and conversion factors for a variable based on its length, scale, and Hubble parameter dependencies.
    
        Args:
            vardescription (str, optional): A description of the variable. Defaults to 'text describing variable'.
            aFact (float, optional): The exponent factor for the expansion scale (a). Defaults to 1.0.
            hFact (float, optional): The exponent factor for the Hubble parameter (h). Defaults to 1.0.
    
        Returns:
            dict: A dictionary containing:
                - 'VarDescription': Description of the variable.
                - 'aexp-scale-exponent': Exponent for the scale factor a.
                - 'h-scale-exponent': Exponent for the Hubble parameter h.
        '''
        return {
            'VarDescription': vardescription, 
            'aexp-scale-exponent': aFact, 
            'h-scale-exponent': hFact
        }
        