fontsize     = 20
# Required libraries
import numpy as np
import re
import mendeleev
import matplotlib.pyplot as plt
from collections import OrderedDict
import h5py 
ElementNames = ['Hydrogen', 'Helium', 'Carbon', 'Oxygen', 'Silicon', 'Iron']
class Elements:
    """
	This class sets the physical parameters for the desired chemical elements that are needed\
        to be able to compute the properties of its absorption lines"
        This class reads properties of transitions from the atom.dat file that comes with VpFit.
        
        Parameters
        ----------    
    
           atomfile  : str
                hdf5 ile that was generated with SpecWizard_atomfile. (default = 'atome_dat.hdf5')
    """
    
    def __init__(self, 
                 atomfile="atom_dat.hdf5"):
        
        # file containing physical parameters of ions and elements
        self.atomdat       = atomfile

        # grep parameters
        self.numeric = '[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
        self.rx      = re.compile(self.numeric, re.VERBOSE)
        
                 
    def Info(self):
        #
        print("This class sets the physical parameters for the desired chemical elements that are needed\
        to be able to compute the properties of its absorption lines")
        
    def ElementParameters(self, ElementNames=ElementNames):
        """ 
	For all desired element, read desired values from atom.dat
         The structure returned is of the form
        (Name, atomic weight, (ionization states))
        Name is the name of the element, e.g. Hydrogen
        atomic weight is the mass of the element in amu (isotopes are not accounted for)
        (ionization stages) is a tuple, containing the considered ionisation states, e.g. (0) for HI, or (0,1) for HeI and HeII
        lines: for each ionization stage, a list of tuples, containing wavelength in Angstrom, and f-value
               these are read from Vpfit's "atom.dat"
        """
        
        elementnames = ElementNames
        parameters   = {}
        # astronomical convention for ions, e.g. HI == neutal hydrogen, HII is ionized hydrogen, etc.
        IonState     = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 
                        'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI', 'XVII',
                       'XVIII', 'XIX', 'XX', 'XXI', 'XXII', 'XXIII', 'XXIV', 'XXV', 'XXVI']
        atomfile = h5py.File(self.atomdat ,'r')
        for element in elementnames:
            try:
                element_ions = atomfile[element].keys()
                element_chsym = atomfile[element].attrs["Symbol"]
                if element_chsym == 'D':
                    weight       = getattr(mendeleev,'H').atomic_weight
                else:
                    weight       = getattr(mendeleev,element_chsym).atomic_weight
                nstates  = len(element_ions)
                states   = {}
                for ion in element_ions:
                    name  = ion
                    states[name] = {}
                    lambda0 = atomfile[element][ion]['lambda0'][...]    
                    fvalue  = atomfile[element][ion]['f-value'][...]
                    states[name]["Lines"] = {'Lambda0': lambda0, 'f-value': fvalue}
            except:
                print('Element not included in the atomfile')

            parameters[element] = {'Weight':weight, 'Nstates':nstates, 'States':states}
        atomfile.close()            
        return parameters                            
 
     
    def Plot(self, ElementNames=ElementNames):
        # plot all transitions used
        elementnames = ElementNames
        parameters   = self.ElementParameters(ElementNames=ElementNames)
        #
        fig, ax  = plt.subplots(1, 1, figsize = (15, 8))
        #
        nelements = len(elementnames)
        plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.CMRmap(np.linspace(0., 1., nelements+1))))
        cycle    = plt.rcParams['axes.prop_cycle'].by_key()['color']
        nsys     = np.arange(nelements)
        
        # wavelength range to plot
        lmin =  100
        lmax = 2400
        for (isys, element) in zip(nsys, elementnames):
            color      = cycle[isys]
            #
            name       = element
            pars       = parameters[element]
            weight     = pars["Weight"]
            ionstates  = pars["States"]
            nionstates = len(ionstates)
            print("Element = {}, weigth = {}, ion stages = {}".format(name, weight, nionstates))
            linestyles = OrderedDict(
                    [('solid',               (0, ())),
                     ('loosely dotted',      (0, (1, 10))),
                     ('dotted',              (0, (1, 5))),
                     ('densely dotted',      (0, (1, 1))),

                     ('loosely dashed',      (0, (5, 10))),
                     ('dashed',              (0, (5, 5))),
                     ('densely dashed',      (0, (5, 1))),

                     ('loosely dashdotted',  (0, (3, 10, 1, 10))),
                     ('dashdotted',          (0, (3, 5, 1, 5))),
                     ('densely dashdotted',  (0, (3, 1, 1, 1))),

                     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
                     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
                     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))])
            
            for ion, (lname, linestyle) in zip(ionstates, enumerate(linestyles.items())):
                #
                try:
                    lambda0 = ionstates[ion]["Lines"]["Lambda0"]  # rest wavelength
                    f_value = ionstates[ion]["Lines"]["f-value"]
                    indxs   = np.arange(len(lambda0))
                    for indx in indxs:
                        w = lambda0[indx]
                        ax.vlines(w, isys, 1+isys, color=color, linestyle=linestyle[1])
                        if indx == 0:
                            label = ion
                            ax.vlines(w, isys, 3+isys, color=color, linestyle=linestyle[1], label=label)
                except:
                    pass
        ax.set_xlim(lmin, lmax)
        ax.set_ylim(0,12)
        ax.legend(frameon=False, fontsize=fontsize, ncol=5)
        ax.set_xlabel(r"$\lambda\ [\AA]$", fontsize=fontsize)
        
        fig.show()
