import unyt
from unyt.dimensions import length, time
from unyt import accepts, returns
import numpy as np
import os

class HIDamping:
    ''' Compute the damping constant for HI Lyman-series  lines '''
    def __init__(self, verbose=False, nmax=5):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        self.EinsteinA_file = os.path.join(script_dir, 'EinsteinA.dat')    # data file with Einstein coefficients
        self.nmax           = nmax              # maximum level to read in
        self.verbose        = verbose
    
        #
        self.A              = self.ReadTabulatedEinsteinCoefficients()

    def ReadTabulatedEinsteinCoefficients(self):
        '''
          Read tabulated Einstein coefficients
          Use these to compute the casecade
        '''

        nmax    = self.nmax
        verbose = self.verbose
        
        # Nist value of forbidden 2s-1s transition. This value is not in the data file read here
        A_2s_1s = 2.496e-06 / unyt.s
        
        
        # columns are n_low, l_low, n_up, l_up, A [1/s]
        import time
        tinit    = time.time()
        dtype    = {'names': ('nd', 'ld', 'nu', 'lu', 'A'),
                  'formats': (np.int32, np.int32, np.int32, np.int32, np.float64)}
        data    = np.loadtxt(self.EinsteinA_file, delimiter=",", dtype=dtype, comments='#', ).T
        tinit   = time.time() - tinit
        if verbose:
            print(" ... Read numerical data in time {0:1.2f}".format(tinit))
    
        # create Einstein coefficients dictionary
        t0      = time.time()
        A       = {}
        # loop over upper level
        for nu in np.arange(2, nmax+1):
            for lu in np.arange(nu):
                conf_i = self.Config(n=nu, l=lu)
                A[conf_i] = {}
                # loop over lower level
                for nd in np.arange(nu):
                    for ld in np.arange(nd):
                        conf_k = self.Config(n=nd, l=ld)
                        A[conf_i][conf_k] = 0
        t0 = time.time() - t0
        if verbose:
            print(" ... Created dictionary of Einstein coefficients in a time {0:1.2f}".format(t0))
                        
        # insert the values from the file
        t1       = time.time()
        nups     = data['nu'][:]
        lups     = data['lu'][:]
        nds      = data['nd'][:]
        lds      = data['ld'][:]
        Avals    = data['A'][:]
        for nup, lup, nd, ld, Aval in zip(nups, lups, nds, lds, Avals):
            conf_i = self.Config(n=nup, l=lup)
            conf_k = self.Config(n=nd, l=ld)
            if nup <= nmax:
                A[conf_i][conf_k] = Aval / unyt.s
            else:
                continue
        t1 = time.time() - t1
        if verbose:
            print(" ... Inserted numerical values in Einstein dictionary in a time {0:1.2}".format(t1))
        
        # insert A_2s-1s
        nu = 2
        lu = 0
        nd = 1
        ld = 0
        conf_i = self.Config(n=nu, l=lu)
        conf_k = self.Config(n=nd, l=ld)
        A[conf_i][conf_k] = A_2s_1s
    
        return A
    
        
    def Config(self, n=1, l=1):
        '''
              configuration states are tuples of the form (n,l), where:
          n = principle quantum number, n=1->nmax
          l = angular momentum number, l=0->n-1
        '''
        return (n,l)
        
    def DeConfig(self, config='1s'):
        '''
            extract n and l value for a given configuration state
        '''
        return config[0], config[1]

    @returns(1/time)
    def DampingConstant(self, n=2):
        '''
        Evaluate and return the damping constant of the n->1 Lyman series transition
        '''
        if self.nmax <= n:
            return 0/ unyt.s
            
        gamma = 0.0
        for nlow in np.arange(1,n):
            for llow in np.arange(0,nlow):
                gamma += self.A[(n, 1)][(nlow, llow)]
        if self.verbose:
            print("n = {}, gamma = {}".format(n, gamma))
        return gamma
        
    @returns(length)    
    def VacuumWavelength(self, n=2):
        '''
        Compute and return the wavelength of the n->1 Lyman series transition
        '''
        #

        wave_n  = 1./(unyt.rydberg_constant * (1. - 1./n**2))
        return wave_n
        
    @returns(length/time)        
    def DampingVelocity(self, n=2):
        '''
        Compute and return the width of the damping wing of the n->1 Lyman series transition in velocity units
        '''
        if self.nmax <= n:
            return 0 * unyt.km / unyt.s
        deltav = self.DampingConstant(n=n) * self.VacuumWavelength(n=n) / (4*np.pi)
        return deltav

