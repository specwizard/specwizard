from SpecWizard_BuildInput import Build_Input
from SpecWizard_Input import ReadData
from SpecWizard_ProjectData import SightLineProjection
from SpecWizard_ComputeOpticaldepth import ComputeOpticaldepth


#class SpecWizard:
    
#    def __init__(self,):
#        a=3
        
        
def GenerateShortSpectra(Wizard=[]):

    snapshot  = ReadData(wizard = Wizard)
    data      = snapshot.read_particles()
    sightlineprojection  = SightLineProjection(Wizard)
    projected_LOS = sightlineprojection.ProjectData(data)

    cspec          = ComputeOpticaldepth(Wizard)
    opticaldepth   = cspec.MakeAllOpticaldepth(projected_LOS)

    return opticaldepth 