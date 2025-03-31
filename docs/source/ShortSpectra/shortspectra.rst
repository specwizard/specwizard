.. _example_page:

Guide for calculating a ShortSpectra
=================
In this guide we will go through the basic process of using `specwizard.GenerateShortSpectra` for the calculation of the optical depth of a simulation sightline.

In order to run this example you would need to have `specwizard` installed, the line transition file and a ionization table. And some simulation data that can be downloaded by:

And the input file that contains the details of what we want to be obtained. Here we have input configuration for a EAGLE full snapshot and an EAGLE line of sight file. 



.. tabs::

    .. tab:: Wizard.yml (LOS) 

        .. code:: YAML

            file_type:
                
                sim_type: eagle
                snap_type: snapshot
                
            snapshot_params:
                
                directory: /cosma7/data/Eagle/ScienceRuns/Planck1/L0100N1504/PE/REFERENCE/data/snapshot_019_z001p004
                file: snap_019_z001p004.0.hdf5

            sightline:
                ProjectionAxes: ['simx','simy','simz']
                ProjectionStart: [0.5,0.5,0]
                ProjectionLength: 1
                SightLength: null
                ProjectionExtend:
                    extend: False
                    extendfactor: 3
                nsight: 0
                
            ionparams:
                table_type: specwizard_cloudy
                iondir: /cosma/home/dp004/dc-aram1/workdir7/pyspecwizard_main/spwtables/HM12/
                fname: null
                ions: [['Hydrogen', 'H I'],['Helium','He II'],['Carbon','C III'],['Carbon','C IV'],['Oxygen','O VI'],['Silicon','Si IV']]
                SFR_properties:
                    modify_particle: True
                    ignore_particle: False
                    Temperature [K]: 1.0E+4    # Putting the plus sign between E and the number is important to not be processed as a string. 
                atomfile: atomic_info.hdf5
                
            ODParams:
                VelOffset_kms: 0
                PecVelEffectsOff: False
                ThermalEffectsOff: False
                VoigtOff: False

    .. tab:: Wizard.yml (snapshot)

        .. code:: YAML

            file_type:
                
                sim_type: eagle
                snap_type: snapshot
                
            snapshot_params:
                
                directory: /cosma7/data/Eagle/ScienceRuns/Planck1/L0100N1504/PE/REFERENCE/data/snapshot_019_z001p004
                file: snap_019_z001p004.0.hdf5

            sightline:
                ProjectionAxes: ['simx','simy','simz']
                ProjectionStart: [0.5,0.5,0]
                ProjectionLength: 1
                SightLength: null
                ProjectionExtend:
                    extend: False
                    extendfactor: 3
                nsight: 0
                
            ionparams:
                table_type: specwizard_cloudy
                iondir: /cosma/home/dp004/dc-aram1/workdir7/pyspecwizard_main/spwtables/HM12/
                fname: null
                ions: [['Hydrogen', 'H I'],['Helium','He II'],['Carbon','C III'],['Carbon','C IV'],['Oxygen','O VI'],['Silicon','Si IV']]
                SFR_properties:
                    modify_particle: True
                    ignore_particle: False
                    Temperature [K]: 1.0E+4    # Putting the plus sign between E and the number is important to not be processed as a string. 
                atomfile: atomic_info.hdf5
                
            ODParams:
                VelOffset_kms: 0
                PecVelEffectsOff: False
                ThermalEffectsOff: False
                VoigtOff: False

        
For more information about the content we refer you to the BuildInput documentation page. 

.. include:: ShortSpectra.ipynb
   :parser: myst_nb.docutils_