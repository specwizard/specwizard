import numpy as np
from SpecWizard_Elements import Elements
from SpecWizard_Lines import Lines
import Phys
constants = Phys.ReadPhys()
#
import scipy.interpolate as interpolate

class ComputeOpticaldepth:
    ''' Methods to compute optical depth as a function of velocity for a single sight line '''
    def __init__(self, sightlineprojection):
        self.specparams = sightlineprojection
        self.periodic   = self.specparams['extra_parameters']['periodic']
        
        # for each of the ions, determine rest-wavelength and f-value of the first transition
        self.elements = self.specparams["elementparams"]["ElementNames"]

        #
        self.transitions = self.specparams["ionparams"]["transitionparams"]
        
        #
        self.constants  = Phys.ReadPhys()
        
        self.ThermEff   = self.specparams['ODParams']['ThermalEffectsOff']
        self.PecVelEff  = self.specparams['ODParams']['PecVelEffectsOff']
        self.VoigtOff   = self.specparams['ODParams']['VoigtOff']
    def MakeAllOpticaldepth(self, projected_los):
        ''' Apply MakeOpticaldepth to compute optical depths for all desired ionic transitions
        Input: 
           projection (dict): output from 
           projection['SightInfo']  = snapshot.ReadParticles(sightline=sightline)["SightInfo"]
           projection['Projection'] = SightLineProjection()

        '''
        
        projection               = {}
        projection["SightInfo"]  = self.specparams["sightline"]
        projection["Header"]     = self.specparams["Header"]
        projection["Projection"] = projected_los
        
        # header information from the snapshot
        self.header = projection["Header"]
        
        # pixel properties
        pixel_kms = self.ToCGS(self.header, projection["Projection"]["pixel_kms"]) / 1e5  # one pixel in km/s
        pixel     = self.ToCGS(self.header, projection["Projection"]["pixel"]) # size of pixel in cm
        sight_kms = self.ToCGS(self.header, projection["SightInfo"]["sightkms"]) / 1e5 
        
            
        # add some extra space of length dv_kms to start and end of spectrum
        # Note: first pixel has vel_kms[0]=0, last pixel has vel_kms[-1]=box_kms-pixel_kms
        vel_kms = np.arange(projection["Projection"]["npix"] ) * pixel_kms
        vunit   = self.SetUnit(vardescription='Hubble velocity',
                                       Lunit=1e5, aFact=0, hFact=0)
        npix    = len(vel_kms)
        #
        
        sightparams = [sight_kms,vel_kms,pixel_kms,pixel] 
        Ions        = self.specparams["ionparams"]["Ions"]
        projectionIW   = projection["Projection"]["Ion-weighted"]
        spectra  = self.WrapSpectra(Ions,projectionIW,sightparams,vel_mod=False,therm_mod=False)
        #return spectra
        if self.ThermEff:

            spectra['ThermalEffectsOff']      = self.WrapSpectra(Ions,projectionIW,sightparams,vel_mod=False,therm_mod=True)

        if self.PecVelEff:
            spectra['PeculiarVelocitiesOff']  = self.WrapSpectra(Ions,projectionIW,sightparams,vel_mod=True,therm_mod=False)
        
        if self.ThermEff and self.PecVelEff:
            spectra['ThermalPecVelOff']       = self.WrapSpectra(Ions,projectionIW,sightparams,vel_mod=True,therm_mod=True)

        DoSimIons   = False
        
        try:
            projectionSIW = projection["Projection"]['SimIon-weighted']
            SimIons       = list(projectionSIW.keys())
            DoSimIons     = True
        except:
            pass
            
        if False:

            SimIons   = np.array(SimIons)
            all_ions      = np.array(Ions)[:,1]
            intsc     = np.in1d(all_ions,SimIons)
            SimElmsIons  = np.array(Ions)[intsc]
            SimElmsIons = [tuple(SimElmsIons[i]) for i in range(len(SimElmsIons))]
            spectra["SimIons"] =  self.WrapSpectra(SimElmsIons,projectionSIW,sightparams,vel_mod=False,therm_mod=False)
            
            if self.ThermEff:
                spectra["SimIons_ThermalEffectsOff"]  = self.WrapSpectra(SimElmsIons,projectionSIW,sightparams,vel_mod=False,therm_mod=True)
        
            if self.PecVelEff:
                spectra['SimIons_PeculiarVelocitiesOff']  = self.WrapSpectra(SimElmsIons,projectionSIW,sightparams,vel_mod=True,therm_mod=False)

            if self.ThermEff and self.PecVelEff:
                spectra['SimIons_ThermalPecVelOff']   = self.WrapSpectra(SimElmsIons,projectionSIW,sightparams,vel_mod=True,therm_mod=True)
            
        return spectra
    
    def WrapSpectra(self,Ions,projection,sightparams,vel_mod=False,therm_mod=False ):
        header     = self.header
        spectra    = {}
        vel_kms    = sightparams[1]
        vunit   = self.SetUnit(vardescription='Hubble velocity',
                               Lunit=1e5, aFact=0, hFact=0)
        for ion in Ions:
            (element_name, ion_name) = ion
            # print("Projecting: ", element_name, ion_name)
            weight  = self.transitions[ion_name]["Mass"] * self.constants["amu"]
            lambda0 = self.transitions[ion_name]["lambda0"]
            f_value = self.transitions[ion_name]["f-value"]
            if lambda0 > 0:
                #
                nions    = self.ToCGS(header, projection[ion_name]["Densities"]) / weight
                vions    = self.ToCGS(header, projection[ion_name]["Velocities"]) / 1e5
                Tions    = self.ToCGS(header, projection[ion_name]["Temperatures"])
                if vel_mod:
                    vions = np.zeros_like(vions)
                if therm_mod:
                    Tions = np.zeros_like(Tions)+1
                
                spectrum = self.MakeOpticaldepth(
                    sightparams=sightparams,
                    weight=weight, lambda0=lambda0, f_value=f_value,
                    nions=nions, vions_kms=vions, Tions=Tions,element_name = element_name)
#                return spectrum
                spectra[ion]                    = {}
                spectra[ion]["Velocities"]      = {'Value': vel_kms, 'Info': vunit}
                spectra[ion]["Optical depths"]  = spectrum["Optical depths"]
                spectra[ion]["Densities"]       = spectrum["Densities"]
                spectra[ion]["Velocities"]      = spectrum["Velocities"]
                spectra[ion]["Temperatures"]    = spectrum["Temperatures"]
                spectra[ion]["TotalIonColumnDensity"] = spectrum["TotalIonColumnDensity"]
                spectra[ion]["Mass"]            = weight
                spectra[ion]["lambda0"]         = lambda0
                spectra[ion]["f-value"]         = f_value
                print("DONE ONE")
        return spectra 
    def MakeOpticaldepth(self, sightparams = [0.0,[0.0],1.0,1.0],
                     weight=1.67382335232e-24, lambda0=1215.67, f_value=0.4164, 
                     nions = [0.0], vions_kms = [0.0], Tions = [0.0],element_name = 'Hydrogen'):

        ''' Compute optical depth for a given transition, given the ionic density, temperature and peculiar velocity '''

        box_kms, vel_kms, pixel_kms,pixel = sightparams
        npix         = len(vel_kms)
        tau          = np.zeros_like(vel_kms)
        densities    = np.zeros_like(vel_kms)
        velocities   = np.zeros_like(vel_kms)  
        temperatures = np.zeros_like(vel_kms)

        
        # passing the extent of the box in km/s introduces periodic boundary conditions
        lines = Lines(v_kms = vel_kms, box_kms=box_kms, constants = self.constants, verbose=False, 
                 lambda0_AA=lambda0, f_value=f_value, naturalwidth_kms=-1,periodic=self.periodic)
            
        # convert from density to column density
        #print("non-zero nions",np.count_nonzero(nions))
        ioncolumns = nions * pixel                                           # in ions/cm^2
        #print("non-zero ioncoluns",np.count_nonzero(nions))

        total_column_density = np.sum(ioncolumns)                            # in ions/cm^2
        dunit        = self.SetUnit(vardescription="Total ion column density", 
                                             Lunit=1.0, aFact=0.0, hFact=0.0)
        total_column_density    = {'Value': total_column_density, "Info": dunit}
#        Tions=np.zeros_like(Tions)+0.001

        # compute b-parameter
        bions_kms = np.sqrt(2*self.constants["kB"]*Tions/weight) / 1e5
        
        # add Hubble velocity to peculiar velocity
        vHubble_kms    = box_kms * np.arange(len(vions_kms)) / len(vions_kms)
        voffset_kms    = self.specparams['ODParams']['Veloffset']  #Default = 0 
        vions_tot_kms  = vions_kms + vHubble_kms + voffset_kms
        

        for (ioncolumn, bion_kms, vion_tot_kms, vion_kms, Tion) in zip(ioncolumns, bions_kms, vions_tot_kms, vions_kms, Tions):
            if ioncolumn > 0:
#                u = lines.EvenBetter(b_kms=bion_kms, v0_kms=vion_tot_kms)
#                return u 
#                print("ioncolumn",ioncolumn)
                dtau          = ioncolumn * lines.sigma * lines.UniversalErf(b_kms=bion_kms, v0_kms=vion_tot_kms)
                tau          += dtau
                densities    += dtau * ioncolumn
                velocities   += dtau * vion_kms
                temperatures += dtau * Tion

        # optical depth-weighted quantities
        if not self.VoigtOff: #element_name == 'Hydrogen':
            tau = lines.ConvolveLorentz(tau)
        mask                = tau > 0
        densities[mask]    /= tau[mask]
        densities          /= pixel
        velocities[mask]   /= tau[mask]
        temperatures[mask] /= tau[mask]

        #
        dunit        = self.SetUnit(vardescription="Tau weighted ion mass density", 
                                             Lunit=1.0, aFact=0.0, hFact=0.0)
        densities   *= weight # mass density
        densities    = {'Value': densities, "Info": dunit}
        vunit        = self.SetUnit(vardescription="Tau weighted ion peculiar velocity",
                                             Lunit=1e5, aFact=0.0, hFact=0.0)
        velocities   = {'Value': velocities, 'Info': vunit}
        tunit        = self.SetUnit(vardescription="Tau weighted ion temperature",
                                             Lunit=1, aFact=0.0, hFact=0.0)
        temperatures = {'Value': temperatures, 'Info': tunit}
        tauunit      = self.SetUnit(vardescription="Ionic optical depth", 
                                             Lunit=1, aFact=0.0, hFact=0.0)
        tau          = {'Value':tau, 'Info':tauunit}
        return {'Optical depths':tau, 'Densities': densities, 'Velocities': velocities, 'Temperatures': temperatures, 'TotalIonColumnDensity':total_column_density}

    def CGSunit(self, header, variable):
        ''' 
        Use the information in the variable to compute the factor needed to convert
        simulation values to proper, h-free cgs units.
        This is of the form
        proper value = simulation value * CGSunit, where
        CGSunit = CGSConversionFactor * h**hexpo * a**aexpo
        '''
        #dependence on expansion factor
        ascale     = (1./(1+header["Cosmo"]["Redshift"]))**variable["Info"]["aexp-scale-exponent"]

        # dependence on hubble parameter
        hscale     = (header["Cosmo"]["HubbleParam"])**variable["Info"]["h-scale-exponent"]
        
        #
        return variable["Info"]["CGSConversionFactor"] * ascale * hscale
        
        
    def ToCGS(self, header, variable):
        ''' 
        return simulations values for this variable in proper cgs units (no h)
        '''
        return variable["Value"] * self.CGSunit(header, variable)

    def SetUnit(self, vardescription = 'text describing variable', Lunit=constants['Mpc'], aFact=1.0, hFact=1.0):
        return {'VarDescription': vardescription, 'CGSConversionFactor':Lunit, 'aexp-scale-exponent' :aFact, 'h-scale-exponent': hFact}

            
            