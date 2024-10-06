
# required python libraries
import numpy as np
import scipy.special as special
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import matplotlib.pyplot as plt
fontsize = 20

''' We define some oft used interpolation kernels '''
class Bspline:
    ''' The B-spline used in Gadget-2 '''
    def __init__(self, nbins=1000):
        self.name  = 'Gadget-2'
        self.norm  = 8 / np.pi
        self.nbins = nbins
        
    def Info(self):
        #
        print("This class defines the B-spline kernel of Monaghan and Lattanzio, as used in Gadget-2")
        
    def Kernel(self,q):
        # Input: q: distance between points in units of smoothing length
        # Output: spline kernel of Monaghan and Lattanzio, as used in Gadget-2
        result = np.zeros_like(q)
        
        # outer boundary
        mask         = (q > 1)
        result[mask] = 0
        
        # outer shell        
        mask         = ((q< 1) & (q>1/2))
        result[mask] = 2*self.norm*(1-q[mask])**3
        
        # inner shell
        mask         = (q<= 1/2)
        result[mask] = self.norm * (1-6*q[mask]**2+6*q[mask]**3)
        return result
    
    def Spline(self):
        # create polynomial interpolation of the 1D spline
        qs     = np.arange(0,1,1./self.nbins)
        ws     = self.Kernel(qs)
        spline = interpolate.interp1d(qs, ws, kind='linear', axis=- 1, copy=True, bounds_error=False, fill_value=0)
        return spline

class TGauss:
    ''' A truncated Gaussian, with dispersion sigma=1. It is 
        truncated at N sigma in *each of the 3 Cartesian directions* '''
    
    def __init__(self, nbins=100):
        self.name       = 'Gaussian'
        self.nsigma     = 3                                                             # extent of kernel in units of sigma
#         self.N          = np.sqrt(2*np.pi) * special.erf(self.nsigma/np.sqrt(2))
#         self.sigmaoverh = np.pi**(1./3.) / (2*self.N)                                   # relation betwen sigma and h
        self.Noverh     = (np.pi)**(1./3.) / 2.                                             # N/h
        self.sigmaoverh = self.Noverh / np.sqrt(2.*np.pi) / special.erf(self.nsigma/np.sqrt(2))  # sigma/N
        self.nbins      = nbins
        
    def Info(self):
        #
        print("This class defines a truncated, normalized Gaussian interpolation kernel")

    def Kernel(self,q):
        # Input: q: distance between points in units of sigma
        result = np.zeros_like(q)
        
        # outer boundary
        mask         = (q > self.nsigma)
        result[mask] = 0
        
        #
        mask         = (q <= self.nsigma)
        result[mask] = np.exp(-q**2/(2*self.sigmaoverh**2)) / (self.Noverh)**3
        return result

    def Spline(self):
        # create polynomial interpolation of the 1D spline
        qs     = np.arange(0, self.nsigma, 1. / self.nbins)
        ws     = self.Kernel(qs)
        spline = interpolate.interp1d(qs, ws, kind='linear', axis=- 1, copy=True, bounds_error=False, fill_value=0)
        return spline
    
    
kernel = Bspline()
kernel.Info()
#kernel = TGauss()

class ColumnTable:
    ''' This class creates a table to interpolate column densities computed through
        an input spherical 3D kernel '''
    def __init__(self, kernel=kernel):
        self.name   = kernel.name
        self.nbins  = kernel.nbins
        self.kernel = kernel.Kernel
        self.spline = kernel.Spline()
        
    def Info(self):
        #
        print("This class uses an input kernel class to create a 2D interpolation table, returning the column density as a function of impact parameter\
        and integration limit. This 2D table can be interpolated, as illustrated in the Plot definition")

    def Column(self):
        ''' Return grid and values for interpolation of the column density of the kernel, 
            the function Column(b, z)'''
        bs      = np.arange(0,1.2,1./self.nbins)   # range of impact parameters
        zs      = np.arange(-1.2,1.2,1./self.nbins)  # range of upper limits of the z-integral
        values  = np.zeros((len(bs), len(zs)))     # will contain the values of column density integration, Column(b, z)
        indxs   = np.arange(len(bs))
        # 
        for (ind,b) in zip(indxs, bs):
            q     = np.sqrt(b**2+zs**2)
            w     = self.kernel(q)
            integ = integrate.cumtrapz(w, zs, initial=0)
            values[ind][:] = integ
        return (bs, zs), values

    def Plot(self):
        ''' test some aspects and illustrates the usage of Kernel and ColumnTable '''
        # returns: a plot
        kernel   = self.kernel       # the kernel, W(q)
        spline   = self.spline       # the spline fit to W(q)

        # generate an instance of the column table
        (btable, ztable), values = self.Column()

        # Verify normalization of the 3D integral
        qs     = np.arange(0,3,1/100)
        func   = 4 * np.pi * qs**2 * kernel(qs)
        integ  = np.trapz(func,qs)
        print('Spherical integral yields {0:1.0f}'.format(integ))

        # plot kernel, and column density as a function of impact parameter
        fig, ax  = plt.subplots(1, 3, figsize = (22, 6))

        #left: spherical kernel, and spline fit
        ax[0].plot(qs, kernel(qs), 'b',label='Kernel')
        qspline = np.arange(0,1.1, 1.1/20)
        kspline = spline(qspline)
        ax[0].plot(qspline, kspline, 'rD', label='spline fit')
        # annotations
        ax[0].set_xlim(0,1.5)
        ax[0].set_ylim(0,2.6)
        ax[0].set_xlabel(r"$q$", fontsize=fontsize)
        ax[0].set_ylabel(r"$W(q,h=1)$", fontsize=fontsize)
#        ax[0].set_title(self.name, fontsize=fontsize)
        ax[0].legend(frameon=False, fontsize=fontsize)

        # right: column-density integral for various values of the impact parameter
        bpars = np.arange(0,1,1./10)   # range of impact parameters
        zvals = np.arange(-2,2,1./20)  # upper limit of column density integration
        for b in bpars:
            bvals  = np.zeros_like(zvals) + b         # impact parameter
            pts    = np.column_stack((bvals, zvals))
            column = interpolate.interpn((btable, ztable), values, pts, bounds_error=False, fill_value=None)
            label  = r"$b={0:1.2f}$".format(b)
            ax[1].plot(zvals, column, label=label)
            diff     = (np.roll(column, -1) - column) / 1
            diff[-1] = 0
            ax[2].plot(zvals,diff,label=label)
        # annotations
        ax[1].set_xlim(-1.5,1.5)
        ax[1].set_ylim(0, 4)
        ax[1].set_xlabel(r"$z$", fontsize=fontsize)
        ax[2].set_xlabel(r"$z$", fontsize=fontsize)

        ax[1].set_ylabel(r"${\rm Column}(b,z)$ (Cumulative)", fontsize=fontsize)
        ax[2].set_ylabel(r"${\rm Column}(b,z)$ (Difference)", fontsize=fontsize)

        ax[1].set_title(self.name, fontsize=fontsize)
        ax[2].legend(frameon=False, fontsize=fontsize-2,handlelength=1)
        
        fig.show()            