
import numpy as np
import Phys
import importlib
import glob
from SpecWizard_Elements import Elements
import h5py
from  reading_simulations import *
Phys = importlib.reload(Phys)
from SimulationInputKeys import get_simkeys

# physical constants in cgs units
constants  = Phys.ReadPhys()


wizard = {}

class ReadData:
    '''
    Read SPH particles from Eagle (Gadget), Hydrangea(C-Eagle) or Swift output
    input: Wizardparamas(output from running SpecWizard_BuildInput)
    '''
    def __init__(self, wizard=wizard):
        self.wizard = wizard
        self.fdir       = wizard["snapshot_params"]["directory"]
        self.fname      = self.fdir + '/' + wizard["snapshot_params"]["file"]
        self.simtype    = wizard["file_type"]["sim_type"]
        self.snaptype   = wizard["file_type"]["snap_type"]
        self.readIonFrac = wizard['extra_parameters']['ReadIonFrac']['ReadIonFrac'] 
        sim_keys  = get_simkeys(self.simtype)
        groupdic  = sim_keys[self.snaptype]
        self.groupdic = groupdic
        groupname = groupdic['groupname']
        self.groupname = groupname    
        self.wizard['short-LOS'] = False
    
    def read_particles(self, ):
        simtype = self.simtype
        if simtype == "eagle":            
            eagle               =  ReadEagle(self.wizard)
            self.header         =  eagle.read_header()
            particles,sightline =  eagle.read_particles()
 
        elif (simtype == "swift") or (simtype == "colibre"):
        
            swift               =  ReadSwift(self.wizard)
            self.header         =  swift.read_header()
            particles,sightline =  swift.read_particles()

        elif simtype == "illustris":
              
            illustris           =  ReadIllustris(self.wizard)
            self.header         =  illustris.read_header()
            particles,sightline =  illustris.read_particles()

        elif simtype == "hydrangea":
              
            hydrangea           =  ReadHydrangea(self.wizard)
            self.header         =  hydrangea.read_header()
            particles,sightline =  hydrangea.read_particles()

        else:
            print("Simulation not supported!") 

            
        inputfunc   = InputFunctions(fileparams=self.wizard,header=self.header)
        self.ToCGS  = inputfunc.ToCGS
        self.Hubble = inputfunc.Hubble
        
        self.wizard['sightline'] = sightline
        self.wizard['Header'] = self.header

        return {'SightInfo': sightline, 'Particles' : particles, 'Header': self.header}