# Required libraries
import importlib
import numpy as np
import re
import mendeleev
from SpecWizard_Elements import Elements
import h5py
import scipy.interpolate as interpolate

elements  = [('Hydrogen',mendeleev.H.atomic_weight,'H0+'), ('Helium',mendeleev.He.atomic_weight,'He1+')]
iondir    = './HM12Tables'

ElementNames = ['Hydrogen', 'Helium', 'Oxygen', 'Carbon', 'Silicon']
atomfile     = "/cosma/home/dphlss/tt/atom.dat"       # ascii file containing physical parameters of transitions, from VpFit
elements     = Elements()

class IonizationBalance:
    def __init__(self, iondir=iondir):
        self.iondir    = iondir  # directory that contains cloudy interpolation tables
        
    def Info(self):
        #
        print("This class computes the fractional abundance of a given ion, by interpolating\
               a table genetrated by Cloudy - assuming photo-ionization equilibrium")

 
    def IonAbundance(self, redshift=0.0, nH_density=1.0, temperature=1e4, ionname = 'H I'):
        ''' Return the fractional abundance of a given ion
              
        Compute and return the ionization levels for the given element
        Inputs:
          redshift     : current redshift
          density:     : proper hydrogen density in particles / cm^3 of the gas
          temperature  : temperature of the gas in K
          ionname      : name of ion,e.g. 'H I' for neutral hydrogen
        Output: log_10 n_ion/n_element
        '''
        
        # OPen the file for interpolation (Create this file using the Cloudy notebook)
        fname = self.iondir + '/' + ionname + '.hdf5'       

        # read the table
        try:
            (table_z, table_LogTs, table_LognHs), table = self.ReadIonTable(fname)
        except:
            print("error: could not read table {} ".format(fname))
            return -1
        
        TInterpol  = np.log10(temperature)
        nHInterpol = np.log10(nH_density)
        zInterpol  = redshift
        pts        = np.column_stack((zInterpol, TInterpol, nHInterpol))
        result     = interpolate.interpn((table_z, table_LogTs, table_LognHs), table, pts, method='linear', bounds_error=False, fill_value=None)
        return result
    
    def ReadIonTable(self, fname="./HM12Tables/H I.hdf5"):
        ''' Read the hdf5 file containing the Cloudy abundances 
         Read table with ionization states for elements
         input: file to read
         output: the tuple (z, LogT, LognH, LogAbundance) where
           z:            tabulated redshifts
           logT:         tabulated values of log temperature in K
           LognH:        tabulated values of log hydrogen number density in units 1/cm^3
           LogAbundance: log aboundance of this ion in units of total, e.g. for ion H+, log n_{HI}/n_H
        '''
        hf  = h5py.File(fname,  "r")
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

    