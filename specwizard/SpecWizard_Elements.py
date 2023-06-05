fontsize     = 20
# Required libraries
import importlib
import numpy as np
import re
import mendeleev
import matplotlib.pyplot as plt
from collections import OrderedDict

ElementNames = ['Hydrogen', 'Helium', 'Carbon', 'Oxygen', 'Silicon', 'Iron']
class Elements:
    ''' 
        This class reads properties of transitions from the atom.dat file that comes with VpFit.
        Input:
           ElementNames: list of elements whose proerties we want to read, for example
               ElementNames = ['Hydrogen', 'Helium']
           atomfile = atome.dat file that comes with VpFit
    '''
    
    def __init__(self, 
                 atomfile="/cosma/home/dphlss/tt/atom.dat"):
        
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
        ''' For all desired element, read desired values from atom.dat
         The structure returned is of the form
        (Name, atomic weight, (ionization states))
        Name is the name of the element, e.g. Hydrogen
        atomic weight is the mass of the element in amu (isotopes are not accounted for)
        (ionization stages) is a tuple, containing the considered ionisation states, e.g. (0) for HI, or (0,1) for HeI and HeII
        lines: for each ionization stage, a list of tuples, containing wavelength in Angstrom, and f-value
               these are read from Vpfit's "atom.dat"'''
        
        elementnames = ElementNames
        Parameters   = {}
        # astronomical convention for ions, e.g. HI == neutal hydrogen, HII is ionized hydrogen, etc.
        IonState     = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 
                        'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI', 'XVII',
                       'XVIII', 'XIX', 'XX', 'XXI', 'XXII', 'XXIII', 'XXIV', 'XXV', 'XXVI']
        
        for element in elementnames:
            if element == 'Hydrogen':
                weight       = mendeleev.H.atomic_weight                
                Nstates      = 2
                States       = {}
                for state in np.arange(Nstates):
                    name                  = 'H {}'.format(IonState[state])  # ion identifier
                    States[name]          = {}
                    States[name]["Lines"] = self.LineSeries(name)

            if element == 'Helium':
                weight  = mendeleev.He.atomic_weight
                Nstates = 3
                States  = {}
                for state in np.arange(Nstates):
                    name         = 'He{}'.format(IonState[state])  # ion identifier
                    States[name] = {}
                    if (state == 0) or (state == 2):
                        States[name]["Lines"] = self.LineSeries(name)
                    if state == 1:
                        # atom.dat does not contain HeII
                        lambda0       = [303.7822]
                        fvalue        = [0.416]
                        States[name]["Lines"] = {'Lambda0': np.array(lambda0), 'f-value': np.array(fvalue)}
                        

            if element == 'Carbon':
                weight  = mendeleev.C.atomic_weight
                Nstates = 6
                States  = {}
                for state in np.arange(Nstates):
                    name                  = 'C {}'.format(IonState[state])    # ion identifier
                    States[name]          = {}
                    States[name]["Lines"] = self.LineSeries(name)
                
            if element == 'Oxygen':
                weight  = mendeleev.O.atomic_weight
                Nstates = 8
                States  = {}
                for state in np.arange(Nstates):
                    name                  = 'O {}'.format(IonState[state])         # ion identifier
                    States[name]          = {}
                    States[name]["Lines"] = self.LineSeries(name)
                
            if element == 'Silicon':
                weight  = mendeleev.Si.atomic_weight
                Nstates = 14
                States  = {}
                for state in np.arange(Nstates):
                    name                  = 'Si{}'.format(IonState[state])         # chemical element name
                    States[name]          = {}
                    States[name]["Lines"] = self.LineSeries(name)
                    
            if element == 'Iron':
                weight  = mendeleev.Fe.atomic_weight
                Nstates = 26
                States  = {}
                for state in np.arange(Nstates):
                    name                  = 'Fe{}'.format(IonState[state])         # chemical element name
                    States[name]          = {}
                    States[name]["Lines"] = self.LineSeries(name)

            if element == 'Magnesium':
                weight  = mendeleev.Mg.atomic_weight
                Nstates = 12
                States = {}
                for state in np.arange(Nstates):
                    name                  = 'Mg{}'.format(IonState[state])         # chemical element name
                    States[name]          = {}
                    States[name]["Lines"] = self.LineSeries(name)

            
            Parameters[element] = {'Weight':weight, 'Nstates':Nstates, 'States':States}
            
        return Parameters

    def LineSeries(self, Identifier):
       # read Lyman series from atom.dat
        fname = self.atomdat
        with open(fname) as fp:
            line     = fp.readline()
            Wave     = []
            Fvalue   = []
            while line:
                line    = fp.readline()
                if re.search(Identifier, line):
                    vals    = self.rx.findall(line)
                    vals    = np.asarray(vals, float)
                    Wave.append(vals[0])
                    Fvalue.append(vals[1])
        Result = {'Lambda0': Wave, 'f-value': Fvalue}
        return Result
    
    def IonValues(self, ion_element=[('Hydrogen',' H I')]):
        ''' extract physical parameters of this ion from the element class '''
        elementpars = self.ElementParameters()
        element,ion = ion_element
        try:
            weight  = elementpars[element]["Weight"]
            lambda0 = elementpars[element]["States"][ion]["Lines"]["Lambda0"]
            fvalue  = elementpars[element]["States"][ion]["Lines"]["f-value"]
        except:
            weight  = -1
            lambda0 = -1
            fvalue  = -1
        return {"Weight": weight, 'Lambda0': lambda0, "f-value":fvalue}


    def HILymanSeries(self):
        # read Lyman series from atom.dat
        with open(fname) as fp:
            line     = fp.readline()
            Wave     = []
            Fvalue   = []
            while line:
                line    = fp.readline()
                if re.search('H I ', line):
                    vals    = self.rx.findall(line)
                    vals    = np.asarray(vals, float)
                    Wave.append(vals[0])
                    Fvalue.append(vals[1])
        Result = {'Lambda0': Wave, 'f-value': Fvalue}

        return Result
    
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
