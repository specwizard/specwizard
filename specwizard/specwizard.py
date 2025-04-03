from .SpecWizard_BuildInput import Build_Input
from .SpecWizard_Input import ReadData
from .SpecWizard_ProjectData import SightLineProjection
from .SpecWizard_ComputeOpticaldepth import ComputeOpticaldepth
from .SpecWizard_Longspectra import LongSpectra
from .SpecWizard_SaveOpticaldepth import OpticalDepth_IO
from .SpecWizard_AnalyseOpticaldepth import Analyse_Opticaldepth


def GenerateShortSpectra(Wizard=[]):

    snapshot  = ReadData(wizard = Wizard)
    data      = snapshot.read_particles()
    sightlineprojection  = SightLineProjection(Wizard)
    projected_LOS = sightlineprojection.ProjectData(data)
    to_physical   = snapshot.to_physical 
    cspec          = ComputeOpticaldepth(Wizard)
    opticaldepth   = cspec.MakeAllOpticaldepth(projected_LOS)
    method          = {}
    method['to_physical'] = to_physical

    opticaldepth['Methods'] = method
    projected_LOS['Methods'] = method 
    data['Methods'] = method
    return opticaldepth,projected_LOS,data