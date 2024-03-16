
import numpy as np
import scipy.interpolate as interpolate
#
from SpecWizard_Input import ReadData
from SpecWizard_Elements import Elements
from SpecWizard_IonizationBalance import IonizationBalance
from SpecWizard_SplineInterpolation import ColumnTable
from SpecWizard_SplineInterpolation import Bspline, TGauss
from SpecWizard_IonTables_test import IonTables
#
import Phys
constants = Phys.ReadPhys()
kernel = Bspline()
kernel = TGauss()

class SightLineProjection:
    
    ''' Interpolate particle properties to a sight line, using SPH kernel interpolation '''
    def __init__(self, specparams, kernelprojection="Bspline",pixkms=1):
        self.specparams        = specparams

        self.kernelprojection = specparams["extra_parameters"]['Kernel']
        self.pixkms           = specparams["extra_parameters"]['pixkms']
        #set kernel
        if kernelprojection=="Bspline":
            columntable      = ColumnTable(Bspline())
            self.kernelprojection = columntable.Column()
        
        elif kernelprojection=="TGauss":
            columntable      = ColumnTable(TGauss())
            self.kernelprojection = columntable.Column()

        else:
            print("ERROR "+kernelprojection+ " is not a valid kernel." + "\n" + "The valid kernels are: "
                  + "\n" + " Bspline" + "\n" + " TGauss")
        
        
        #set periodic
        self.periodic          = self.specparams['extra_parameters']['periodic']
        
        #set pixel size
        self.specparams["pixkms"] = self.pixkms
    def Info(self):
        #
        print("This class interpolates particle data to a sightline")
    
    def PeriodicDist(self, dx, period, periodic=True):
        # impose peridic boundary conditions, with period = 1
        if periodic:
            dx[dx<0.5 * period] += period
            dx[dx>0.5 * period] -= period
        return dx
    
    def PixelSizeCheck(self,pixel_size,hmin):
        '''This compares the smoothing length and the pixel size to check if 
           you are under/over shoothing with the resolution
           Input: Pixel Size, Minimum smoothing length
           Output: A Warning.        
        '''
        ratio = pixel_size  / hmin
        
        if ratio >= 2.0:
            print("Warning! the pixel size is {:0.1f} times larger than the smallest smoothing lenght. Probably you are not resolving nicely high density peaks".format(1/ratio))
        elif 1/ratio <= 2.0:
            print("Warning! the pixel size is {:0.1f} times smaller than the smallest smoothing lenght.Probably you are overshooting with the pixel size".format(ratio))

    def CalculateHDI(self,header,Hfrac,Mgas):
        ''' We use the procedure showed in F. van de Voort et all 2012 (eq 2) for the calculation of deuterium abundance (D/H)
            
        
        '''

        Minitial      = self.ToCGS(header,header['MassTable']['GasMass'])
        Mgas          = self.ToCGS(header,Mgas)
        Primordial_H  = 0.76
        Primordial_DH = 2.547E-5 
        
        DAbund        = (Minitial/Mgas) * (Primordial_H/Hfrac) * Primordial_DH
        
        return DAbund
    
    def ProjectData(self, sightlinedata,ReadIonfrac=False):
        ''' Interpolate particle properties to a sight line'''
        sightinfo          = sightlinedata["SightInfo"]
        particles          = sightlinedata["Particles"]
        ions               = self.specparams["ionparams"]["Ions"]
        ionizationbalance  = self.specparams["ionparams"]["IonizationBalance"]
        elementnames       = self.specparams["elementparams"]["ElementNames"]
        header             = sightlinedata["Header"] 
        los_length         = sightinfo['ProjectionLength']
        
        
        boxkms = sightinfo["Boxkms"]["Value"]  # in km/s
        box    = sightinfo["Box"]["Value"]     # in programme units
        self.specparams['short-LOS']     = sightinfo['short-LOS']
        
        # compute extra properties for sightline
        
        # ensure that an integer number of pixels fit into a sight line
        sightkms = boxkms * los_length
        npix    = np.array((sightkms / self.specparams["pixkms"]),dtype=int) + 1
        pixkms  = sightkms / npix
        sight   = box * los_length
        pix     = sight / npix
        
        # Sight line properties
        # number of pixels of sight line and z-values of pixels
        zpix  = np.arange(-npix, 2*npix) / (3*npix.astype(float))  # extend the sightline for easy periodic boundary implementation
        
        # (x,y) positions of sight line as a fraction of the box 
        proj0 = xproj = sightinfo['x-position'] * sightinfo["Box"]["Value"]
        proj1 = yproj = sightinfo['y-position'] * sightinfo["Box"]["Value"]
        proj2 = zproj = sightinfo['z-position'] * sightinfo["Box"]["Value"]

        
        # +++++ Particle properties
        # mass
        mass  = particles['Masses']['Value']
        
        # smoothing length and its inverse
        h     = particles['SmoothingLengths']['Value']
        hinv  = 1./h
        
        # Check smoothing length vs pix size 
        #self.PixelSizeCheck(pix,h.min())
        

        # impact parameter in units of smoothing length. Note: we impose periodic boundary conditions if periodic=True.
        dx  = self.PeriodicDist(particles['Positions']['Value'][:,0] - proj0, box,self.periodic)*hinv  # off-set from sightline in x
        dy  = self.PeriodicDist(particles['Positions']['Value'][:,1] - proj1, box,self.periodic)*hinv  # off-set from sightline in x
        b   = np.sqrt(dx**2+dy**2)                                                         # impact parameter
        vz  = particles['Velocities']['Value'][:,2]# peculiar velocity along sightline
            
       # mask particles that contribute
        mask  = b < 1.1

        #shift_in_z is 
        shift_in_z   = zproj 
        int_zmins    = ((particles['Positions']['Value'][:,2] - h -shift_in_z) / pix).astype(int) - 1
        zcents       = particles['Positions']['Value'][:,2] - shift_in_z
        int_zcents   = np.round(zcents / pix).astype(int)
        int_zmaxs    = ((particles['Positions']['Value'][:,2] + h - shift_in_z) / pix).astype(int) + 1


        
        # interpolated values along the sight line
        
        # total density
        rho_tot          = {}
        rho_tot['Densities']     = {'Value': np.zeros(npix), 'Info': particles['Densities']['Info']}      # density
        rho_tot['Velocities']    = {'Value': np.zeros(npix), 'Info': particles['Velocities']['Info']}     # density-weighted peculiar velocity
        rho_tot['Temperatures']  = {'Value': np.zeros(npix), 'Info': particles['Temperatures']['Info']}   # density-weigthed temperatue

        # element densities
        rho_element             = {}
        nunit                   = particles['Densities']['Info']
        nunit['VarDescription'] = 'Element mass-densities'
        vunit                   = particles['Velocities']['Info']
        vunit["VarDescription"] = 'Element-weighted velocities'
        tunit                   = particles['Temperatures']['Info']
        tunit["VarDesciption"]  = 'Element-weighted temperatures '

        # properties per element
        for element in elementnames:
            rho_element[element]                 = {}
            rho_element[element]['Densities']    = {'Value': np.zeros(npix), 'Info': nunit}   # element mass density                          # ion particle density
            rho_element[element]['Velocities']   = {'Value': np.zeros(npix), 'Info': vunit}   # element-weighted peculiar velocity
            rho_element[element]['Temperatures'] = {'Value': np.zeros(npix), 'Info': tunit}   # element-weighted temperature
            rho_element[element]['Mass']         = self.specparams["elementparams"][element]["Mass"] 
        
        # ion densities
        rho_ion                 = {}
        nunit                   = particles['Densities']['Info']
        nunit['VarDescription'] = 'Ion mass-densities'
        vunit                   = particles['Velocities']['Info']
        vunit["VarDescription"] = 'Ion-weighted velocities'
        tunit                   = particles['Temperatures']['Info']
        tunit["VarDesciption"]  = 'Ion-weighted temperatures '

        # variables per ion
        for (element, ion) in ions:
            rho_ion[ion]                 = {}
            rho_ion[ion]['Densities']    = {'Value': np.zeros(npix), 'Info': nunit}  # ion mass density
            rho_ion[ion]['Velocities']   = {'Value': np.zeros(npix), 'Info': vunit}  # ion-weighted peculiar velocity
            rho_ion[ion]['Temperatures'] = {'Value': np.zeros(npix), 'Info': tunit}  # ion-weigthed temperature
            rho_ion[ion]['Mass']         = self.specparams["ionparams"]["transitionparams"][ion]["Mass"]
            rho_ion[ion]['lambda0']      = self.specparams["ionparams"]["transitionparams"][ion]["lambda0"]
            rho_ion[ion]['f-value']      = self.specparams["ionparams"]["transitionparams"][ion]["f-value"]

            
        # determine element fractions
#        if ReadIonfrac==False:
        ParticleAbundances = {}
        hydrogenfraction = self.ToCGS(header, particles["Abundances"]["Hydrogen"])                                                                  
        nH_cgs           = self.ToCGS(header, particles['Densities']) * hydrogenfraction / constants["mH"]
        temperature      = self.ToCGS(header, particles['Temperatures'])
        Z                = self.ToCGS(header, particles['Metallicities'])
        redshift         = header["Cosmo"]["Redshift"] + np.zeros_like(temperature)
        for element in elementnames:
            massfraction   = self.ToCGS(header, particles["Abundances"][element])
            ParticleAbundances[element] = {}
            ParticleAbundances[element]["massfraction"] = massfraction  # fraction of this element by mass
            ParticleAbundances[element]["element mass"] = \
                self.specparams["elementparams"][element]["Mass"] * constants["amu"] # mass of element [g]

        # determine ion fractions
        ComputedIonFractions = {}
        for (element, ion) in ions:
            
            if ion == 'D I':
                ComputedIonFractions[ion] = self.CalculateHDI(header, hydrogenfraction, particles["Masses"])# n_ion/n_element
            else:
                LogIonfraction = ionizationbalance.IonAbundance(redshift, 
                                    nH_density=nH_cgs, temperature=temperature,Z=Z, ion=ion)
                ComputedIonFractions[ion] = 10.**LogIonfraction # n_ion/n_element
        
        try:
            SimulationIonFraction = {}
            SimIons = np.array(list(particles['SimulationIonFractions'].keys()))
            rho_ion_sim           = {}
            for SimIon in SimIons:

                SimulationIonFraction[SimIon] = self.ToCGS(header, particles['SimulationIonFractions'][SimIon])

                rho_ion_sim[SimIon]                 = {}
                rho_ion_sim[SimIon]['Densities']    = {'Value': np.zeros(npix), 'Info': nunit}  # ion mass density
                rho_ion_sim[SimIon]['Velocities']   = {'Value': np.zeros(npix), 'Info': vunit}  # ion-weighted peculiar velocity
                rho_ion_sim[SimIon]['Temperatures'] = {'Value': np.zeros(npix), 'Info': tunit}  # ion-weigthed temperature
                rho_ion_sim[SimIon]['Mass']         = self.specparams["ionparams"]["transitionparams"][SimIon]["Mass"]
                rho_ion_sim[SimIon]['lambda0']      = self.specparams["ionparams"]["transitionparams"][SimIon]["lambda0"]
                rho_ion_sim[SimIon]['f-value']      = self.specparams["ionparams"]["transitionparams"][SimIon]["f-value"]
            
        except:
            pass
        # interpolation table for projected kernel
        (btable, ztable), table = self.kernelprojection  

       # particle loop
        npart = len(h)
        for i in np.arange(npart):
            if b[i] > 1: # require particle contributes to sightline
                continue
            
            zcent   = zcents[i]                                                  # z-location of particle
            izmin   = int_zmins[i]                                               # z-coordinate of first pixel that particle contributes to (maybe negative)
            izmax   = int_zmaxs[i]                                               # z-coordinate of first pixel that particle contributes to (maybe larger than sight line)

            # print(zcent, self.lospart["pos"][i,2])
            zvals    = (np.arange(izmin, izmax).astype(float)) * pix              # z-values of pixels that this particle contributes to
            intz     = np.arange(izmin, izmax)                                    # corresponding pixel coordinates
 
            zvals   -= zcent                                                      # z-distance between pixel and particle
            zvals   *= hinv[i]                                                    # z-distance beween pixel and particle in units of smoothing length
            
            # restrict to pixels that are inside the smoothing length
            mask     = (zvals**2 + b[i]**2) < 1.1
            zvals    = zvals[mask]                                          
            intz     = np.array(intz[mask])
            

            if len(zvals) == 0 or len(zvals) == 1:
                # print("Error len zvals = 0 .")
                # print("Pix = ", pix)
                # print("smoothing length= ", h[i])
                # print(izmin,izmax)
                # print("zpos",zcent)
                # print("assigned pix", int_zcents[i] * pix)
#                zvals = np.array(int_zcents[i] * pix) + h[i]
                pts   = np.array([b[i],ztable.max()])
                diff  = interpolate.interpn((btable, ztable), table, pts, bounds_error=True, fill_value=None)
                diff *= hinv[i]**2
                diff /= pix 
                diff *= mass[i]
                intz  = int_zcents[i]
            
            else:
                bvals    = np.zeros_like(zvals) + b[i]                                # impact parameter in units of smoothing length
                pts      = np.column_stack((bvals, zvals))                            # points (b,z) that particle contributes to
    #            column   = interpolate.interpn((btable, ztable), table, pts, bounds_error=False, fill_value=None) # cumulative column density contributed to these pixels
                try:
                    column   = interpolate.interpn((btable, ztable), table, pts, bounds_error=True, fill_value=None)
                except:
                    print("bounds error: ", izmin, izmax, 1./hinv[i])
                    print("btable: ", btable.min(), btable.max())
                    print("ztable: ", ztable.min(), ztable.max())
                    print("bvals = ", bvals)
                    print("zvals = ", zvals)
                    
                column  *= hinv[i]**2 # scale integral for h=1, to the actual value of h
                
                
                # difference the cumulative column to get the contribution to each pixel
                diff     = (np.roll(column, -1) - column) / pix
                diff[-1] = 0
                
                # account for periodic boundary conditions
                if self.periodic:
                    intz[intz<0]    += npix
                    intz[intz>=npix]-= npix
                else:
                    mask = (intz >=0) & (intz < npix)
                    intz = intz[mask]
                    diff = diff[mask]
                    
                diff  *= mass[i]       # multiply with mass of particle                                                 
            
            mdiff  = np.copy(diff) # save mass-weighted kernel contribution
            
            # density-weighted quantities
            rho_tot['Densities']['Value'][intz]     += diff
            rho_tot['Velocities']['Value'][intz]    += diff * vz[i]
            rho_tot['Temperatures']['Value'][intz]  += diff * temperature[i]

            # Densities of elements
            for element in elementnames:
                diff = mdiff * ParticleAbundances[element]["massfraction"][i]
                #
                rho_element[element]['Densities']['Value'][intz]    += diff
                rho_element[element]['Velocities']['Value'][intz]   += diff * vz[i]
                rho_element[element]['Temperatures']['Value'][intz] += diff * temperature[i]
        
            # Ion densities
            for (element, ion) in ions:
                diff = mdiff * ParticleAbundances[element]["massfraction"][i] * ComputedIonFractions[ion][i]
                #
                rho_ion[ion]['Densities']['Value'][intz]    += diff
                rho_ion[ion]['Velocities']['Value'][intz]   += diff * vz[i]
                rho_ion[ion]['Temperatures']['Value'][intz] += diff * temperature[i]
            
            try:
                Ions2do = np.array([ions[i][1] for i in range(len(ions))])
                maskindx = np.where(np.in1d(Ions2do,SimIons))[0]

                for (element, ion) in np.array(ions)[maskindx]:
                    diff = mdiff * ParticleAbundances[element]["massfraction"][i] * SimulationIonFraction[ion][i]
                    rho_ion_sim[ion]['Densities']['Value'][intz]    += diff
                    rho_ion_sim[ion]['Velocities']['Value'][intz]   += diff * vz[i]
                    rho_ion_sim[ion]['Temperatures']['Value'][intz] += diff * temperature[i]

            except:
                pass
         
        # normalize peculiar velocity and temperature
        mask = rho_tot['Densities']['Value'] > 0
        rho_tot['Velocities']['Value'][mask]   /= rho_tot['Densities']['Value'][mask]
        rho_tot['Temperatures']['Value'][mask] /= rho_tot['Densities']['Value'][mask]
        
        for element in elementnames:
            mask = rho_element[element]['Densities']['Value'] > 0
            rho_element[element]['Velocities']['Value'][mask]    /= rho_element[element]['Densities']['Value'][mask]
            rho_element[element]['Temperatures']['Value'][mask]  /= rho_element[element]['Densities']['Value'][mask]

        for (element, ion) in ions:            
            mask = rho_ion[ion]['Densities']['Value'] > 0
            rho_ion[ion]['Velocities']['Value'][mask]    /= rho_ion[ion]['Densities']['Value'][mask]
            rho_ion[ion]['Temperatures']['Value'][mask]  /= rho_ion[ion]['Densities']['Value'][mask]

        
        # prepare output
        unit                   = particles["Positions"]['Info']
        unit["VarDescription"] = 'pixel size'
        pixelsize              = {'Value': pix, 'Info': unit}
        pixelsize_dv           = self.SetUnit(vardescription='Hubble velocity accross pixel', 
                                              Lunit=1e5,
                                              aFact=0.0,
                                              hFact=0.0)
        pixelsize_kms          = {'Value': pixkms, 'Info': pixelsize_dv}
        
        boxkms_info            = sightinfo["Boxkms"]['Info']
        sightkms               = {'Value':sightkms,'Info':boxkms_info}

        self.specparams['sightline']['sightkms'] = sightkms
        
        result = {'pixel': pixelsize, 'pixel_kms': pixelsize_kms,
                  'npix': npix,
                  'Header': header, 
                  'Element-weighted': rho_element,
                  'Ion-weighted': rho_ion,
                  'Mass-weighted': rho_tot}
                  

        try:
            
            for SimIon in SimIons:
                mask = rho_ion_sim[SimIon]['Densities']['Value'] > 0
                rho_ion_sim[SimIon]['Velocities']['Value'][mask] /= rho_ion_sim[SimIon]['Densities']['Value'][mask]
                rho_ion_sim[SimIon]['Temperatures']['Value'][mask]  /= rho_ion_sim[SimIon]['Densities']['Value'][mask]
            
            result['SimIon-weighted'] = rho_ion_sim
        except:
            pass

        
        return result
    
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