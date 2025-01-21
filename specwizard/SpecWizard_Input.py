import numpy as np
import glob
import h5py
from .Phys import ReadPhys
from .SpecWizard_Elements import Elements
from .reading_simulations import *
from .SimulationInputKeys import get_simkeys

# Physical constants in CGS units
constants = ReadPhys()

wizard = {}

class ReadData:
    """
    This class reads SPH particle data from simulation outputs such as Eagle (Gadget),
    Hydrangea (C-Eagle), or Swift. It requires `Wizard` dictionary, generated from running
    `SpecWizard_BuildInput`, to configure its input parameters.

    Args:
        wizard (dict): A dictionary containing all required simulation and file 
            parameters for data reading, snapshot configurations, and additional settings. 

    Attributes:
        wizard (dict): Configuration parameters for the simulation.
        fdir (str): Directory of the snapshot files.
        fname (str): Full path to the snapshot file.
        simtype (str): Type of simulation (e.g., "eagle", "swift").
        snaptype (str): Snapshot type as specified in wizard.
        readIonFrac (bool): Whether to read ion fraction information.
        groupdic (dict): Dictionary of group configurations for the snapshot type.
        groupname (str): Group name for accessing data in the simulation file.
        header (dict): Header information from the simulation snapshot.
        ToCGS (function): Conversion function for CGS units.
        Hubble (function): Hubble parameter based on input configuration.
    """
    
    def __init__(self, wizard=wizard):
        """
        Initialize the ReadData class with a specified `wizard` dictionary.

        Args:
            wizard (dict): Dictionary of simulation parameters for reading data. 
        """
        self.wizard = wizard
        self.fdir = wizard["snapshot_params"]["directory"]
        self.fname = self.fdir + '/' + wizard["snapshot_params"]["file"]
        self.simtype = wizard["file_type"]["sim_type"]
        self.snaptype = wizard["file_type"]["snap_type"]
        self.readIonFrac = wizard['extra_parameters']['ReadIonFrac']['ReadIonFrac']
        sim_keys = get_simkeys(self.simtype)
        groupdic = sim_keys[self.snaptype]
        self.groupdic = groupdic
        groupname = groupdic['groupname']
        self.groupname = groupname    
        self.wizard['short-LOS'] = False
        self.header = self.read_header()
        inputfunc = InputFunctions(fileparams=self.wizard, header=self.header)
        self.ToCGS = inputfunc.ToCGS
        self.Hubble = inputfunc.Hubble
        self.to_unyt = inputfunc.set_particle_data_to_unyt
        self.to_physical = inputfunc.to_physical


    def read_header(self):
        """
        Read and return header information from the simulation snapshot based on its type.

        Returns:
            dict: A dictionary containing header information of the simulation snapshot.

        Raises:
            ValueError: If the simulation type is not supported.
        """
        simtype = self.simtype
        if simtype == "eagle":            
            eagle = ReadEagle(self.wizard)
            header = eagle.read_header()
 
        elif simtype in ["swift", "colibre"]:
            swift = ReadSwift(self.wizard)
            header = swift.read_header()

        elif simtype == "illustris":
            illustris = ReadIllustris(self.wizard)
            header = illustris.read_header()

        elif simtype == "hydrangea":
            hydrangea = ReadHydrangea(self.wizard)
            header = hydrangea.read_header()

        else:
            raise ValueError("Simulation not supported!")
        
        return header

    def read_particles(self):
        """
        Read particle data from the simulation snapshot and return sightline, particle, 
        and header information.

        Returns:
            dict: A dictionary with keys:
                - 'SightInfo' (dict): Information about the sightline.
                - 'Particles' (dict): Particle data from the simulation.
                - 'Header' (dict): Header information from the simulation snapshot.

        Raises:
            ValueError: If the simulation type is not supported.
        """
        simtype = self.simtype
        if simtype == "eagle":            
            eagle = ReadEagle(self.wizard)
            particles, sightline = eagle.read_particles()
 
        elif simtype in ["swift", "colibre"]:
            swift = ReadSwift(self.wizard)
            particles, sightline = swift.read_particles()

        elif simtype == "illustris":
            illustris = ReadIllustris(self.wizard)
            particles, sightline = illustris.read_particles()

        elif simtype == "hydrangea":
            hydrangea = ReadHydrangea(self.wizard)
            particles, sightline = hydrangea.read_particles()

        else:
            raise ValueError("Simulation not supported!")

        self.wizard['sightline'] = sightline
        self.wizard['Header'] = self.header
        particles = self.to_unyt(part_data=particles)

        return {'SightInfo': sightline, 'Particles': particles, 'Header': self.header}
