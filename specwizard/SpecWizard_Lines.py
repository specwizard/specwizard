
import numpy as np
import Phys
constants = Phys.ReadPhys()
#
from scipy.signal import convolve
from scipy.special import erf
from astropy.modeling.functional_models import Voigt1D
from scipy.special import voigt_profile as VoigtSciPy
import scipy.interpolate as interpolate

# convolve(in1, in2, mode='full', method='auto')[source]
class Lines:
    ''' Methods to compute optical depth as a function of velocity for a single absorber,
        implementing:
        - a Gaussian profile
        - a Voigt profile
        '''
    def __init__(self, v_kms =0.0, box_kms=-1.0, constants=constants, verbose=False, 
                 lambda0_AA=1215.67, f_value=0.4164, naturalwidth_kms=6.06076e-3,periodic=True):
        self.constants    = constants
        self.verbose      = verbose
        self.v_kms        = np.array(v_kms)    # velocity bins to evaluate optical depth [km/s]
        self.box_kms      = box_kms            # if >0, periodic extent of the spectrum
        self.box2_kms     = 0.5 * box_kms      # half the box size
        self.lambda0      = lambda0_AA * 1e-8  # rest-wavelength           [cm]
        self.f_value      = f_value            # oscillator strength       [dimensionless]
        self.naturalwidth = naturalwidth_kms   # natural line width        [km/s]
        self.pix_kms      = v_kms[1] - v_kms[0] # pixel size
        # self.naturalwidth = 1.0
        # line cross section times speed of light, c [cm^2 * cm/s]
        self.sigma        = constants["c"] * np.sqrt(3*np.pi*constants["sigmaT"]/8.) * self.f_value * self.lambda0
        self.periodic     = periodic
        # compute cumulative integral of Gaussian with sigma=1.0
        # compute a spline iinterpolation function for this integral
        v_max = 30
        v_min = 0
        self.nfine      = 1000000
        self.du         = (v_max-v_min)/self.nfine
        ukms            = np.arange(0,self.nfine)*self.du
        center          = 0.5*(v_max+v_min) 
        self.CG         = 0.5 * (1.+erf((ukms-center)))
        
    def UniversalErf_old(self,b_kms=10,v0_kms=0.0 ):
        
        
        def FromFineToCoarse(FVal,b,dv):
            nfine  = self.nfine
            du     = self.du
            nfine2 = int(nfine*0.5)
            onestep = int(((FVal-nfine2)*b)+nfine2)
            return int((onestep) * (du/dv))

        def FromCoarseToFine(Cval,b,dv):
            nfine  = self.nfine
            du     = self.du
            finb_eqv = int((dv/du) * Cval)#-displacement
            nfine2 = int(nfine*0.5)
            return int(nfine2+(finb_eqv-nfine2)*1/b)
        auxvk0 = v0_kms
#        try:
        vkms        = self.v_kms
        v0_kms     -= 0.5*vkms.max()
        n_coarse    = len(vkms)
        convVelfact = (30/vkms.max()) #30 is the v_max of the fine universal erf
        bfine       = b_kms * convVelfact
        v0_kms     *= convVelfact
        dv          = (30/n_coarse)
        vkms        = np.arange(0,n_coarse)* dv
        firstOne    = 694314
        lastZero    = 302613
        CoarseFirstOne = FromFineToCoarse(firstOne,bfine,dv)
        CoarseLastZero = FromFineToCoarse(lastZero,bfine,dv)
        CoarseRange    = np.arange(CoarseLastZero,CoarseFirstOne+1)    

        finalCG        = np.zeros_like(vkms)
        vF             = np.vectorize(FromCoarseToFine)
        FineRange      = vF(CoarseRange,bfine,dv)

        finalCG[CoarseRange] = self.CG[FineRange]
        finalCG[:CoarseRange[0]] = 0
        finalCG[CoarseRange[-1]:] = 1
        v0_displacement = (1/dv)*v0_kms
        if v0_displacement < 0:
            v0_displacement = round((1/dv)*v0_kms)
        else:
            v0_displacement = int((1/dv)*v0_kms)

        cphi   = finalCG
        phi    = (np.roll(cphi, -1) - cphi) 
        phi[phi==-1] = 0
        phi     = np.roll(phi,1)
        phi = np.roll(phi,v0_displacement)

        if self.periodic == False:
            print("esta pasando")
            indx_disp = (np.arange(n_coarse) + v0_displacement)                
            phi[indx_disp[indx_disp<0]+n_coarse] = 0
            phi[indx_disp[indx_disp>=n_coarse]-n_coarse] = 0           
#        except:
#            u      = self.PeriodicVelocity(v0_kms=auxvk0+0.5*self.pix_kms) / b_kms
#            phi    = np.exp(-u**2) / (np.sqrt(np.pi) * b_kms) # Gaussian line-profile function
            
        phi /= self.pix_kms
        return phi/1e5

    def UniversalErf(self,b_kms=10,v0_kms=0.0 ):
        
        
        def FromFineToCoarse(FVal,b,dv):
            nfine  = self.nfine
            du     = self.du
            nfine2 = int(nfine*0.5)
            onestep = int(((FVal-nfine2)*b)+nfine2)
            return int((onestep) * (du/dv))

        def FromCoarseToFine(Cval,b,dv):
            nfine  = self.nfine
            du     = self.du
            finb_eqv = int((dv/du) * Cval)#-displacement
            nfine2 = int(nfine*0.5)
            return int(nfine2+(finb_eqv-nfine2)*1/b)
        
        auxvk0 = v0_kms
        try:
            vkms        = self.v_kms
            v0_kms     -= 0.5*vkms.max()
            n_coarse    = len(vkms)
            convVelfact = (30/vkms.max()) #30 is the v_max of the fine universal erf
            bfine       = b_kms * convVelfact
            v0_kms     *= convVelfact
            dv          = (30/n_coarse)
            vkms        = np.arange(0,n_coarse)* dv
            firstOne    = 694314
            lastZero    = 302613
            CoarseFirstOne = FromFineToCoarse(firstOne,bfine,dv)
            CoarseLastZero = FromFineToCoarse(lastZero,bfine,dv)
            CoarseRange    = np.arange(CoarseLastZero,CoarseFirstOne+1)    

            finalCG        = np.zeros_like(vkms)
            if len(CoarseRange)>2:
                vF             = np.vectorize(FromCoarseToFine)
                FineRange      = vF(CoarseRange,bfine,dv)
                finalCG[CoarseRange] = self.CG[FineRange]
                finalCG[:CoarseRange[0]] = 0
                finalCG[CoarseRange[-1]:] = 1
                v0_displacement = (1/dv)*v0_kms
                if v0_displacement < 0:
                    v0_displacement = round((1/dv)*v0_kms)
                else:
                    v0_displacement = int((1/dv)*v0_kms)

                cphi   = finalCG
                phi    = (np.roll(cphi, -1) - cphi) 
                phi[phi==-1] = 0
                phi     = np.roll(phi,1)
                phi = np.roll(phi,v0_displacement)
                phi /= self.pix_kms


            else:
                FineRange   = np.array([len(self.CG)*0.5],dtype=int)
                CoarseRange = np.array([n_coarse*0.5],dtype=int)

                finalCG[CoarseRange] = self.CG[FineRange]
                finalCG[:CoarseRange[0]] = 0
                finalCG[CoarseRange[-1]:] = 1

                cphi   = finalCG
                phi    = (np.roll(cphi, -1) - cphi) 
                phi[phi==-1] = 0
                phi     = np.roll(phi,1)
                phi /= self.pix_kms
                peak_locat  = np.argmax(phi)
                disp  = (np.round(auxvk0 / self.pix_kms) - peak_locat).astype(int)
                phi  = np.roll(phi,disp)

            if self.periodic == False:
                #print("esta pasando")
                indx_disp = (np.arange(n_coarse) + v0_displacement)                
                phi[indx_disp[indx_disp<0]+n_coarse] = 0
                phi[indx_disp[indx_disp>=n_coarse]-n_coarse] = 0           
        except:
            #print("gaussian")
            u      = self.PeriodicVelocity(v0_kms=auxvk0+0.5*self.pix_kms) / b_kms
            phi    = np.exp(-u**2) / (np.sqrt(np.pi) * b_kms) # Gaussian line-profile function
            phi /= self.pix_kms

        return phi/1e5

    
    
    def PeriodicVelocity(self, v0_kms=0.0):
        '''
           Compute and return periodic velocity off-set from line centre

        '''
        # off-set from line centre
        voffset_kms = self.v_kms - v0_kms
        if self.periodic:
        # apply periodic boundary conditions if box_kms>0
            if self.box_kms > 0:
                voffset_kms[voffset_kms >  self.box2_kms] -= self.box_kms
                voffset_kms[voffset_kms < -self.box2_kms] += self.box_kms

        return voffset_kms
    


    
    def DirectGauss(self, b_kms=10.0, v0_kms=0.0):
        ''' return line-profile function for Gaussian profile, Equation 3
        input:
             b_kms:    float: b-parameter [km/s]
             v0_kms:   float: velocity of line centre [km/s]
        output: 
             line profile    : float or array [s/cm]
        '''
        u      = self.PeriodicVelocity(v0_kms=v0_kms+0.5*self.pix_kms) / b_kms
        phi    = np.exp(-u**2) / (np.sqrt(np.pi) * b_kms) # Gaussian line-profile function
        return phi/1e5 # convert to units of [s/cm]
    
    def Lorentz(self, v0_kms=0.0):
        ''' 
        return line-profile function for Lorentian profile, Equation 4
        input:
             v0_kms (float): velocity at line centre [km/s]
        output: 
            line profile (float or array): line shape with unit columnn density [s/cm])
            '''
        width        = self.naturalwidth
        u            = self.PeriodicVelocity(v0_kms=v0_kms) / width
        phi          = (1./(width*np.pi)) / (u**2 + 1.0)               # Lorentzian line-profile function
        return phi/1e5 # convert to units of [s/cm]
        
        
    
    def SciPyVoigt(self, b_kms=10., v0_kms=0.0):
        ''' 
        return Voigt line-profile function, Equation 5
        this version uses the SciPy implementation of the Voigt function
        input:
             b_kms (float):     b-parameter [km/s]
             v0_kms (float):    velocity at line centre [km/s]
        output: 
             line profile (float or array) : line shap with unit column density [s/cm]
        
        '''
        u           = self.PeriodicVelocity(v0_kms=v0_kms) / b_kms
        # 
        sigma_G     = 1.0 / np.sqrt(2.0)                # variance of Gaussian
        gamma_L     = self.naturalwidth / b_kms         # Half-width half-maximum of Lorenzian - the parameter "a"
        # evaluate Voigt profile
        vnorm       = VoigtSciPy(u, sigma_G, gamma_L, out=None)
        phi         = vnorm / b_kms                     # SciPy returns a normalize Voigt profile, which includes the 1/sqrt(pi)
        return phi/1e5  # convert to units of [s/cm]
        
    def Voigt(self, b_kms=10., v0_kms=0.0):
        ''' 
        return Voigt line-profile function, Equation 5
        wrapper of astropy Voigt profile function
        input:
             b_kms (float):  b-parameter [km/s]
             v0_kms (float): velocity at line centre [km/s]
        output: 
            tau (float or array):    optical depth at v_kms [-]
        
        '''
        u           = self.PeriodicVelocity(v0_kms=v0_kms) / b_kms
        fwhm_L      = 2 * self.naturalwidth / b_kms     # FWHM of Lorentzian
        fwhm_G      = 2 * np.sqrt(np.log(2.))           # FWHM of Gaussian
        amplitude_L = 2 / (np.pi * fwhm_L)              # generate unit Voigt
        
        # generate Voigt profile function
        v1    = Voigt1D(x_0=0, amplitude_L= amplitude_L, fwhm_L=fwhm_L, fwhm_G=fwhm_G)
        vnorm = v1(u)
        phi   = vnorm / b_kms # Voigt returns a normalize Voigt profile, which includes the 1/sqrt(pi)
        return phi/1e5 # convert to units of [s/cm]
    
    def ConvolveLorentz(self, phi):
        ''' return convolution of input line profile with Lorentzian
        input: phi(v_z): line profile function as a function of velocity, v_z
        output: phi(v_z) inpout line profile now convolved with natural line profile
        '''
        # The core of the Lorentzian needs to be sampled to 0.01 km/s to get an accurate
        # convolution. Therefore we interpolate the original profile to a finer velocity grid
        # before performing the convolution.
        
        # Create velocity bins for interpolation
        dv         = 1e-3                   # pixel size in km/s
        vmin       = np.min(self.v_kms)
        vmax       = np.max(self.v_kms)
        nbins      = np.int((vmax-vmin)/dv)
        dv         = (vmax-vmin)/float(nbins)
        v_convolve = vmin + np.arange(nbins) * dv 
        
        # Create lorentz profile
        phi_fine    = np.interp(v_convolve, self.v_kms, phi)
        width       = self.naturalwidth
        v_bins      = v_convolve - np.mean(v_convolve)              # centre Lorenz at the centre of the velocity interval
        lorentz     = (1./np.pi) * width / (v_bins**2 + width**2)   
        lorentz     = lorentz / 1e5                                 # convert to units of [s/cm]
        #
        phi_fine    = convolve(phi_fine, lorentz, mode='same') / np.sum(lorentz)
        
        #
        result      = np.interp(self.v_kms, v_convolve, phi_fine)
        
        return result
    
    