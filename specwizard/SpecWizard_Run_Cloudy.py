import mendeleev
import time
import numpy as np
import os
import re
import h5py
import scipy.interpolate as interpolate
import scipy.integrate as integrate
import roman
from .SpecWizard_Elements import Elements
from .Phys import ReadPhys


class Run_cloudy:
	def __init__(self,cloudy_yml='cloudy.yml'):

		with open(cloudy_yml) as file:
			cloudy_yml = yaml.load(file, Loader=yaml.FullLoader)

		cloudycode    = cloudy_yml['cloudycode']

		# Adapt the following line to point to a directory that cloudy can use to store intermediate results. These may be large,
		# so you may want to use /tmp or similar
		cloudyindir   = cloudy_yml['cloudyindir']
		cloudyUVB     = cloudy_yml['cloudyUVB']                                                          # Builtin Haardt & Madau 2012 UV background
		cloudycompo   = cloudy_yml['cloudycompo']                                 # Builtin abundance of elements, no grains
		cloudyscale   = cloudy_yml['cloudyscale']                                                      # 0.01 % of the MW ISM abundance
		cloudyhdf5dir = cloudy_yml['cloudyhdf5dir']                                                    # directory containing .hdf5 tables
		zgrid     = cloudy_yml['zgrid']                                 # redshift range: start, end, step
		LogTgrid  = cloudy_yml['LogTgrid']                          # log10(T[k}]) range: start, end, step
		LognHgrid = cloudy_yml['LognHgrid']                                 # log10(nH(1/cm^3)) total hydrogen number density: start, end, step
		cloudytestgrid   = {'z':zgrid, 'LogT':LogTgrid, 'LognH':LognHgrid}              # cloudygrid of values of (z, LogT and LognH to be computed)
		self.cloudyrun    = {'grid':cloudytestgrid, 'code':cloudycode, 'indir':cloudyindir, 'compo':cloudycompo, 'element scale factor':cloudyscale, 'UVB':cloudyUVB, 'hdf5dir':cloudyhdf5dir} 
		self.elements_to_do = cloudy_yml['elements_to_calculate']
		self.physconst  = ReadPhys()
		self.elements = Elements(cloudy_yml['atom_file'])

	def create_files(self):
		z_range = self.CloudyRange("z")
		for z in z_range:
			name = self.write_cloudy_z_file(z)
			self.run_cloudy(name,z)

	def run_cloudy(self,name,z):
		cloudyrun = self.cloudyrun
		fname = name
		command = 'cd {0:s};export OMP_NUM_THREADS=8; {1:s} -r {2:s}; cd'.format(cloudyrun["indir"],cloudyrun["code"],fname)
		start_time = time.time()
		stream     = os.popen(command)
		print(command)
		output     = stream.read()
		end_time   = time.time()
		duration   = end_time - start_time
		print(" ... Cloudy run finished for z= {0:2.2f}; calculation time = {1:4.1f} seconds".format(z, duration))


	def write_cloudy_z_file(self, z):
		
		Ignore = ['Beryllium', 'Boron', 'Scandinium', 'Lithium', 'Vanadium', 'Copper', 'Fluorine', 
				'Zinc', 'Cobalt', 'Titanium', 'Potassium',         
				'Chlorine', 'Manganese', 'Phosphorus', 'Chromium', 'Nickel', 'Sodium', 
				'Argon']
		#  Save level population for these ions
		# Species = ['H[1:2]', 'He[1:3]', 'C+3[1:3]', 'O+5[1:4]']
		Species = []  # don't specify level population of individual ions

		elements = self.elements

		cloudyrun = self.cloudyrun

		name  = "{0:s}_z{1:2.2f}".format(cloudyrun["UVB"], z)
		fname = '{0:s}'.format(cloudyrun["indir"]) + '/' + name + '.in'
		# write the file
		file  = open(fname, 'w+')
		file.write("table {0:s} redshift {1:2.2f} \n".format(cloudyrun["UVB"], z))
		fname = name + '.cont'
		file.write('save continuum "{0:s}" ## save continuum \n'.format(fname))
		# selected abundance pattern
		file.write('abundances {}\n'.format(cloudyrun["compo"]))
		# elements to be ignored
		for el in Ignore:
			file.write('element off {0:s}\n'.format(el))
		# scale metal abundance by scale factor
		file.write('metals {} log ## reduce abundances of all metals \n'.format(np.log10(cloudyrun["element scale factor"])))
		# Temperature grid
		file.write('constant temperature 4 vary ##do not solve for T equilibrium\n')
		file.write("grid {}  {}  {} ## vary log temperature [K] (start, end, step) \n"
				.format(cloudyrun["grid"]["LogT"]["min"], cloudyrun["grid"]["LogT"]["max"], cloudyrun["grid"]["LogT"]["step"]))
		# density grid
		file.write('hden -5 vary        ## log_10 of H density in cm^-3\n')
		file.write("grid {}  {}  {} ## vary log total H density [particles cm^3] (start, end, step) \n"
				.format(cloudyrun["grid"]["LognH"]["min"], cloudyrun["grid"]["LognH"]["max"], cloudyrun["grid"]["LognH"]["step"]))
		# geometry
		file.write('set dr 0             ## set zone thickness of 1 cm\n')
		file.write('stop zone 1          ## do only one zone\n')
		file.write('iterate to convergence ## iterate\n')
		# files to save
		fname = name + '.heat'
		file.write('save heating last "{0:s}" no hash ## save heating\n'.format(fname))
		fname = name + '.cool'
		file.write('save cooling last "{0:s}" no hash ## save cooling\n'.format(fname))
		fname = name + '.ovr'
		file.write('save overview last "{0:s}" no hash ## save overview\n'.format(fname))
		fname = name + '.ion'
		file.write('save ionization means last "{0:s}" no hash ##save ionization levels\n'.format(fname))
			# save photo-ionization rate for each ion
			#    elementpars = Elements.ElementParameters()
		for element in self.elements_to_do:
				# number of ionization states
			print(element)
			if element=='Sulphur':
				nstates = elements.ElementParameters(['Sulfur'])['Sulfur']["Nstates"]
			else:
				nstates = elements.ElementParameters([element])[element]["Nstates"]
			for state in np.arange(nstates):
				ion   = "{0:1d}".format(state+1)
				fname = name + '_' + element + '_' + ion + '.ionrate'
				# file.write('save ionization rates {} "{}" last\n'.format(element, fname))
				file.write('save gammas element {} {} "{}" \n'.format(element, ion, fname))
			if len(Species) > 1:
				file.write('save species densities last "{0:s}" no hash ## save population levels\n'.format(fname))
				for species in Species:
					file.write('"{0:s}"\n'.format(species))
				#file.write('end\n')
		file.close()

		return name

	def test_cloudy(self):
		command = 'echo "test" | {}'.format(self.cloudyrun["code"])
		stream = os.popen(command)
		output = stream.read()
		print(output)


	def CloudyRange(self,variable='LogT'):
		cloudygrid = self.cloudyrun["grid"]
		minval   = cloudygrid[variable]["min"]
		maxval   = cloudygrid[variable]["max"]
		step     = cloudygrid[variable]["step"]
		result   = np.arange(minval, maxval+step, step)
		return result[result <= maxval+1e-1]

	def write_cloudy_tables_to_hdf5(self,saveheat=True,savecontiniuum=True,saveabundances=True):
		UVB             = self.cloudyrun["UVB"]
		words           = UVB.split(maxsplit=1)
		nwords          = len(words)
		line            = []
		numeric         = '[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
		self.numeric    = numeric
		self.rx         = re.compile(numeric, re.VERBOSE)

		if nwords > 1:
			name = words[0]
			for word in words[1:]:
				name = name + ' ' + word
		else:
			name = UVB
		self.UVB_name = name

		if saveheat:
			self.save_heat()
		if savecontiniuum:
			self.save_continiuum()

		if saveabundances:
			self.save_abundances()
		print("Finished writing hdf5")

		# expression to extract numerical values from strings using regular expression matching

	def save_heat(self):
		CloudyRange = self.CloudyRange
		cloudyrun   = self.cloudyrun
		zs     = np.array(CloudyRange("z"))
		izs    = np.arange(len(zs))
		Htot   = np.zeros( (len(CloudyRange("z")), len(CloudyRange("LogT")), len(CloudyRange("LognH"))))
		Ctot   = np.zeros( (len(CloudyRange("z")), len(CloudyRange("LogT")), len(CloudyRange("LognH"))))
		H1     = np.zeros( (len(CloudyRange("z")), len(CloudyRange("LogT")), len(CloudyRange("LognH"))))
		He1    = np.zeros( (len(CloudyRange("z")), len(CloudyRange("LogT")), len(CloudyRange("LognH"))))
		He2    = np.zeros( (len(CloudyRange("z")), len(CloudyRange("LogT")), len(CloudyRange("LognH"))))
		Comp   = np.zeros( (len(CloudyRange("z")), len(CloudyRange("LogT")), len(CloudyRange("LognH"))))    
		#
		for (iz, z) in zip(izs, CloudyRange("z")):
			heat   = self.ReadHeating(z)
			#
			Htot[iz,:,:] = heat["Htot"][:,:]
			Ctot[iz,:,:] = heat["Ctot"][:,:]
			H1[iz,:,:]   = heat["H1"][:,:]
			He1[iz,:,:]  = heat["He1"][:,:]
			He2[iz,:,:]  = heat["He2"][:,:]
			Comp[iz,:,:] = heat["Comp"][:,:]

		# save all input spectra in hdf5 file
		hdf5name = '{}'.format(cloudyrun["hdf5dir"]) + '/' + 'Heating.hdf5'
		hf       = h5py.File(hdf5name, 'w')
		#store information
		info = hf.create_group("CloudyInfo")
		info.attrs["UVB"] = cloudyrun["UVB"]
		info.attrs["Composition"] = cloudyrun["compo"]
		info.attrs["Element scale factor"] = cloudyrun["element scale factor"]
		# write datasets
		dset = hf.create_dataset('LogT', data=CloudyRange("LogT"))
		dset.attrs['Info'] = "Log_10 of temperature [K]"
		dset = hf.create_dataset('LognH', data=CloudyRange("LognH"))
		dset.attrs['Info'] = "Log_10 of hydrogen number density [1/cm^3]"            
		dset = hf.create_dataset('z', data=CloudyRange("z"))
		dset.attrs['Info'] = "Redshifts [-]"
		#
		dset = hf.create_dataset('Htot', data=Htot)
		dset.attrs['Info'] = "Total heating rate in erg/cm^3/s"
		dset = hf.create_dataset('Ctot', data=Ctot)
		dset.attrs['Info'] = "Total cooling rate in erg/cm^3/s"
		#
		dset = hf.create_dataset('H1', data=H1)
		dset.attrs['Info'] = "Fraction of heating due to H 1"
		dset = hf.create_dataset('He1', data=He1)
		dset.attrs['Info'] = "Fraction of heating due to He 1"
		dset = hf.create_dataset('He2', data=He2)
		dset.attrs['Info'] = "Fraction of heating due to He 2"
		dset = hf.create_dataset('Comp', data=Comp)
		dset.attrs['Info'] = "Fraction of heating due to Comp"    

		#
		hf.close()
		print(" ... written file {}".format(hdf5name))  

	def save_continiuum(self,):
		zs     = np.array(self.CloudyRange("z"))
		izs    = np.arange(len(zs))
		for (iz, z) in zip(izs, zs):
			cont   = self.Continuum(z)
			hnu    = cont["hnu_erg"]
			nufnu  = cont["nufnu_erg/cm2/s"]
			if iz == 0:
				InputSpec = np.zeros((len(self.CloudyRange("z")), len(nufnu)))
			InputSpec[iz,:] = nufnu
		# save all input spectra in hdf5 file
		hdf5name = '{}'.format(self.cloudyrun["hdf5dir"]) + '/' + 'IncidentSpectrum.hdf5'
		hf       = h5py.File(hdf5name, 'w')
		#store information
		info = hf.create_group("CloudyInfo")
		info.attrs["UVB"] = self.cloudyrun["UVB"]
		info.attrs["Composition"] = self.cloudyrun["compo"]
		info.attrs["Element scale factor"] = self.cloudyrun["element scale factor"]
		# write datasets
		dset = hf.create_dataset('z', data=zs)
		dset.attrs['Info'] = "Redshifts [-]"
		dset = hf.create_dataset('hnu', data=hnu)
		dset.attrs['Info'] = "Photon energy hnu [erg]"
		dset = hf.create_dataset('nufnu', data=InputSpec)
		dset.attrs['Info'] = "Incident flux in erg/cm^2/s"
		#
		hf.close()
		print(" ... written file {}".format(hdf5name)) 


	def save_abundances(self):
	# read and store the abundances of each element for each redshift
		elements = self.elements

		CloudyRange = self.CloudyRange
		cloudyrun   = self.cloudyrun
		for element in self.elements_to_do:
			pars      = elements.ElementParameters([element])[element]
			name    =   list(pars["States"].keys())[0].split(" ")[0]
			nlevels   = pars['Nstates']
			result    = np.zeros((len(CloudyRange("z")), len(CloudyRange("LogT")), len(CloudyRange("LognH"))))
			resultne  = np.zeros((len(CloudyRange("z")), len(CloudyRange("LogT")), len(CloudyRange("LognH"))))
			ionrate   = np.zeros(len(CloudyRange("z")))
			for level in np.arange(nlevels):
				file_name = name+" "+roman.toRoman(int(level)+1)
				hdf5name = '{}'.format(cloudyrun["hdf5dir"]) + '/' + file_name + '.hdf5'
				# read output for each redshift
				zs     = np.array(CloudyRange("z"))
				izs    = np.arange(len(zs))
				for (iz, z) in zip(izs, CloudyRange("z")):

					# read the output
					IonLevel   = self.OneIonLevel(z,element, level, nlevels)
					LogTgrid   = IonLevel["LogT"]
					LognHgrid  = IonLevel["LognH"]
					Logne      = IonLevel["Logne"]
					LogLevel   = IonLevel["LogAbundance"]
					# (LogTgrid, LognHgrid, LogLevel) = IonLevel
					result[iz,:,:]   = LogLevel[:,:]
					resultne[iz,:,:] = Logne[:,:]
					ionrate[iz]      = IonLevel["Gamma"] 

				if element == 'Hydrogen':
					# write electron abundance file    
					fname = '{}'.format(cloudyrun["hdf5dir"]) + '/' + 'Electron' + '.hdf5'
					hf    = h5py.File(fname, 'w')
					#store information
					info = hf.create_group("CloudyInfo")
					info.attrs["UVB"] = cloudyrun["UVB"]
					info.attrs["Composition"] = cloudyrun["compo"]
					info.attrs["Element scale factor"] = cloudyrun["element scale factor"]
					# write datasets
					dset = hf.create_dataset('LogT', data=LogTgrid)
					dset.attrs['Info'] = "Log_10 of temperature [K]"
					dset = hf.create_dataset('LognH', data=LognHgrid)
					dset.attrs['Info'] = "Log_10 of hydrogen number density [1/cm^3]"            
					dset = hf.create_dataset('z', data=zs)
					dset.attrs['Info'] = "Redshifts [-]"
					dset = hf.create_dataset('Logne', data=resultne)
					dset.attrs['Info'] = "Log_10 of electron number density [1/cm^3]"
					#
					hf.close()

				# write hdf5 file for this ion
				hf  = h5py.File(hdf5name, 'w')
				#store information
				info = hf.create_group("CloudyInfo")
				info.attrs["UVB"] = cloudyrun["UVB"]
				info.attrs["Composition"] = cloudyrun["compo"]
				info.attrs["Element scale factor"] = cloudyrun["element scale factor"]
				# write datasets
				dset = hf.create_dataset('LogT', data=LogTgrid)
				dset.attrs['Info'] = "Log_10 of temperature [K]"
				dset = hf.create_dataset('LognH', data=LognHgrid)
				dset.attrs['Info'] = "Log_10 of hydrogen number density [1/cm^3]"            
				dset = hf.create_dataset('z', data=zs)
				dset.attrs['Info'] = "Redshifts [-]"
				dset = hf.create_dataset('Gamma', data=ionrate)
				dset.attrs['Info'] = "Photo-ionization rate [1/s]"
				dset = hf.create_dataset('LogAbundance', data=result)
				dset.attrs['Info'] = "Log_10 of fraction of this ionization stage"
				#
				hf.close()
				print(" ... written file {}".format(hdf5name))

	def write_info(self):
		# write some info of each file
		for file in self.files:
			ftype = file[0]
			fname = file[1]
			print("Contents of file {} with extension {} ".format(fname, ftype))
			try:
				with open(fname) as fp:
					line  = fp.readline()
					print('Variables {}'.format(line))
			except:
				pass
	def Continuum(self,z):
		# Parse the cloudy .cont file
		# returns: incident spectrum, in the form of
		#  hnu: photon energy in erg
		#  nufnu: continuum in erg/cm^2/s

		name            = self.UVB_name + "_z{0:2.2f}".format(z)
		basename        = cloudyrun["indir"] + '/' + name
		files = [['ovr', basename + '.ovr'], ['ion', basename + '.ion'], ['cont', basename + ".cont"], ['heat', basename + ".heat"]] # read overview file (.ovr) and ionization stages file (.ion)

		for file in files:
			if (file[0] == 'cont'):
				fname = file[1]
		with open(fname) as fp:
			line      = fp.readline()
			Vars      = line.split()    # list of variable names (e.g. Te, the (electron) temperature)
			nVars     = len(Vars)
			nu        = []
			nufnu     = []
			line      = fp.readline()
			Vars      = line.split()    # list of variable names (e.g. Te, the (electron) temperature)
			nVars     = len(Vars)
			Vals      = []
			gridfound = False
			nlines    = 0
			while line and not gridfound:
				line  = fp.readline()
				vals  = self.rx.findall(line)
				Vars  = line.split() 
				try:
					if Vars[0] == '###########################':
						gridfound = True
				except:
					continue
				else:
					if not gridfound:
						vals  = np.asarray(vals, float)
						nu.append(vals[0])
						nufnu.append(vals[1])
						nlines += 1
		fp.close()

		# save interpolated continuum, converting photon energy to ergs
		hnusave    = np.arange(-1,3,0.01)               # photon energy in Rydberg - saved from 0.1 to 100 Rydberg
		f          = interpolate.interp1d(np.log10(nu), np.log10(nufnu), fill_value="extrapolate", bounds_error=False)   
		nufnusave  = 10**f(hnusave)
		hnusave    = self.physconst['Ryd'] * 10**hnusave

		#
		return {'hnu_erg':hnusave, 'nufnu_erg/cm2/s':nufnusave}
		
	def sigmaHI(self, hnu):
		# HI photo-ionization cross section
		# input: photon energy hnu in erg
		# output: photo-ionization cross section in cm^2/Hz
		# Fit from Verner+'92'
		hnu1    = self.physconst['Ryd']
		aH      = np.sqrt(hnu/hnu1-1)
		result  = 6.3e-18 * (aH**2+1)**(-4) * np.exp(4-4*(np.arctan(aH)/aH)) / (1.-np.exp(-2.*np.pi/aH))    
		if len(hnu) > 1:
			result[hnu<hnu1]= 0
		else:
			if hnu<hnu1:
				result = 0
		return result

	def NameValue(self, line=[], element = 'He', level = 2):
		# Helper routine for ReadHeating
		Vars        = np.array(line.split())
		elementindx = np.where(Vars==element)[0]
		value       = -1 # initialize
		for i in np.arange(len(elementindx)):
			indx = elementindx[i]
			if level >= 0:
				nlev = int(Vars[indx+1])
				if nlev == level:
					value = Vars[indx+2]
			else:
				value = Vars[indx+1]
		return value

	def ReadHeating(self,z):
		# Parse the cloudy .heat file
		# returns: 
		# electron temperature, Te
		# total heating rate, Htot [erg/cm^3/s]
		# total cooling rate, Ctot [erg/cm^3/s]
		# fractional contributions to Htot from H1, He1 and He2 photo-ionization, and from Compton heating (Comp)
		name            = self.UVB_name + "_z{0:2.2f}".format(z)
		basename        = cloudyrun["indir"] + '/' + name
		files = [['ovr', basename + '.ovr'], ['ion', basename + '.ion'], ['cont', basename + ".cont"], ['heat', basename + ".heat"]] # read overview file (.ovr) and ionization stages file (.ion)
		


		for file in files:
			if (file[0] == 'heat'):
				fname = file[1]
		#
		cloudygrid = self.cloudyrun["grid"]        
		LogTs      = self.CloudyRange("LogT")
		LognHs     = self.CloudyRange("LognH")
		nLogTs     = len(LogTs)
		nLognHs    = len(LognHs)
		#
		indx = 0
		Te   = np.zeros((nLogTs, nLognHs))
		Htot = np.zeros((nLogTs, nLognHs))
		Ctot = np.zeros((nLogTs, nLognHs))
		H1   = np.zeros((nLogTs, nLognHs))
		He1  = np.zeros((nLogTs, nLognHs))
		He2  = np.zeros((nLogTs, nLognHs))
		Comp = np.zeros((nLogTs, nLognHs))
		with open(fname) as fp:
			# header line
			line     = fp.readline()
			#
			while line:
				line    = fp.readline()
				if line:
					nHindx  = np.mod(indx, nLognHs)
					Tindx   = ((indx - nHindx) / nLognHs).astype(int)
					Vars    = np.array(line.split())
					Te[Tindx, nHindx]      = float(Vars[1])
					Htot[Tindx, nHindx]    = float(Vars[2])
					Ctot[Tindx, nHindx]    = float(Vars[3])
					H1[Tindx, nHindx]      = self.NameValue(line, element = 'H', level = 1)
					He1[Tindx, nHindx]     = self.NameValue(line, element = 'He', level = 1)
					He2[Tindx, nHindx]     = self.NameValue(line, element = 'He', level = 2)
					Comp[Tindx, nHindx]    = self.NameValue(line, element = 'Comp', level = -1)
					#
					indx += 1
		fp.close()
		return {'Te':Te, 'Htot':Htot, 'Ctot':Ctot, 'H1':H1, 'He1':He1, 'He2':He2, 'Comp':Comp}


	def Overview(self, verbose=False):
		elements = self.elements
		# Parse the cloudy overview output files, making sure we are reading the expected values of T and nH
		LogTs    = self.CloudyRange("LogT")
		LognHs   = self.CloudyRange("LognH")
		nLogTs   = len(LogTs)
		nLognHs  = len(LognHs)
		grid     = np.zeros((nLogTs, nLognHs, 2))
		
		for file in self.files:
			if (file[0] == 'ovr'):
				fname = file[1]

		# make a list of all the variable names in the output file, and a list of all the corresponding values
		with open(fname) as fp:
			line     = fp.readline()
			Vars     = line.split()    # list of variable names (e.g. Te, the (electron) temperature)
			nVars    = len(Vars)
			Vals = []
			while line:
				line  = fp.readline()
				vals  = self.rx.findall(line)
				vals  = np.asarray(vals, float)
				if len(vals >= nVars):
					Vals.append(vals[0:nVars])    # corresponding values
		fp.close()
						
		# Parse variable names to look for index of Te (temperature) and hden (nH)
		indxs = np.arange(len(Vars))
		for (indx, var) in zip(indxs, Vars):
			if var == "Te":
				iT = np.copy(indx)
			if var == "hden":
				inH = np.copy(indx)
			if var == "eden":
				ineden = np.copy(indx)
		print(" T = {}, hden = {}, eden = {}".format(iT, inH, ineden))
		
		# Find indices of Te, and hden (Hydrogen number density, 1/cm^3)
		indxs = np.arange(len(Vals))
		for (indx, Val) in zip(indxs, Vals):
			T        = Val[iT]
			nH       = Val[inH]
			nHindx   = np.mod(indx, nLognHs)
			Tindx    = ((indx - nHindx) / nLognHs).astype(int)

		# Extract range of temperature and densities computed
		LognHVals = Vals[0:nLognHs][nHindx]
		print("LognHvals = {}".format(LognHVals))
			
		# Collect outputs
		Outputs = {}
		Ivars   = np.arange(nVars)
		for (var, ivar) in zip(Vars, Ivars):
			Values = []
			for val in Vals:
				Values.append(val[ivar])
			if verbose:
				print("Collected variable {}".format(var))
			result = {'Variable':var, 'Values':Values}
			Outputs[var] = {'Values':np.array(Values)}
		
		return Outputs


	def OneIonLevel(self, z,element='Hydrogen', level=0, nlevels=3):
		# Parse the .ion file, returning the ionization state of level 'level' for this element
		# Convention is such that level uses the chemical notation
		# Example:
		# for element = 'Hydrogen' and level=0, the function return the 2D grid of log_10 (n_HI/n_H)     .. ie log of the HI fraction
		# for element = 'Carbon', and level=3,  the function return the 2D grid of log_10 (n_C3+/n_C)    .. ie log of the CarbonIV fraction
		#   note: for Hydrogen, level=3 refers to molecular Hydrogen
		name            = self.UVB_name + "_z{0:2.2f}".format(z)
		basename        = cloudyrun["indir"] + '/' + name
		files = [['ovr', basename + '.ovr'], ['ion', basename + '.ion'], ['cont', basename + ".cont"], ['heat', basename + ".heat"]] # read overview file (.ovr) and ionization stages file (.ion)

		PhotoRates = True
		if level >= nlevels-1:
			PhotoRates = False
		Gammaprt = -1
		if PhotoRates:
			fname = basename + '_' + element + '_' + "{}".format(level+1) + '.ionrate'
			with open(fname) as fp:
				line     = fp.readline()
				indx     = 0
				while line and indx < 5:
					line     = fp.readline()                
					Vars     = line.split()    # list of variable names (e.g. Te, the (electron) temperature)
					Vars     = np.array(Vars)
					# identify label "total=""
					ind      = 0
					indtotal = -1
					value    = -1
					for var in Vars:
						if var == 'total=':
							indtotal = ind
							Gammaprt = float(Vars[indtotal+1])  # photo-ionization rate [1/s]
						ind += 1
					indx += 1
			fp.close()

		# Expected grid of temperatures and densities
		cloudygrid = self.cloudyrun["grid"]
		LogTs      = self.CloudyRange("LogT")
		LognHs     = self.CloudyRange("LognH")
		nLogTs     = len(LogTs)
		nLognHs    = len(LognHs)        
		
		for file in files:
			if (file[0] == 'ion'):
				ionfile = file[1]
			if  (file[0] == 'ovr'):
				ovrfile = file[1]
		 
		# parse .ovr file
		fovr  = open(ovrfile)
		line  = fovr.readline()   # read header line
		Vars  = line.split()      # list of variable names (e.g. Te, the (electron) temperature) 
		# identify indx of columns that contain Te and hden
		indxs = np.arange(len(Vars))
		for (indx, var) in zip(indxs, Vars):
			if var == 'Te':
				indxTe = np.copy(indx)
			if var == 'hden':
				indxnH = np.copy(indx)
			if var == 'eden':
				indxeden = np.copy(indx)
				
		# collect values of LogT and LognH 
		LogTvals  = np.zeros((nLogTs, nLognHs))
		LognHvals = np.zeros((nLogTs, nLognHs))
		Lognevals = np.zeros((nLogTs, nLognHs))
		indx      = 0
		line      = fovr.readline()        
		while line:
			vals     = self.rx.findall(line)
			vals     = np.asarray(vals, float)
			LogT     = np.log10(vals[indxTe])
			LognH    = np.log10(vals[indxnH])
			Logeden  = np.log10(vals[indxeden])
			nHindx   = np.mod(indx, nLognHs)
			Tindx    = ((indx - nHindx) / nLognHs).astype(int)
			indx     += 1
			line      = fovr.readline()              
			#
			LogTvals[Tindx, nHindx]  = LogT
			LognHvals[Tindx, nHindx] = LognH
			Lognevals[Tindx, nHindx] = Logeden
		fovr.close()
		LogTgrid  = LogTvals[:,0]
		LognHgrid = LognHvals[0,:]
			
		# parse .ion file
		LogIonvals = np.zeros((nLogTs, nLognHs))
		indx       = 0
		with open(ionfile) as fion:
			line  = fion.readline()
			while line:
				line  = fion.readline()
				words = line.split(maxsplit=1)
				# check if we have the right element
				if len(words) > 1:
					if(words[0] == element):
						numbers      = self.rx.findall(words[1])
						nC           = len(numbers)
						if nC > nlevels:
							nC = nlevels
						levelC       = np.asarray(numbers, float)  # cloudy does not return level if too low
						#
						levels       = -30. + np.zeros(nlevels, dtype=float)
						levels[0:nC] = levelC[0:nC]
						#
						nHindx = np.mod(indx, nLognHs)
						Tindx  = ((indx - nHindx) / nLognHs).astype(int)
						indx  += 1
						#
						LogIonvals[Tindx, nHindx] = levels[level]
		fion.close()
		#
		return {'LogT': LogTgrid, 'LognH':LognHgrid, 'LogAbundance':LogIonvals, 'Logne':Lognevals, 'Gamma':Gammaprt}
