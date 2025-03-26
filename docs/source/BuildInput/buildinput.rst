.. _example_page:

Guide for BuildInput
=================
In order to use ``specwizard`` we need to indicate the program the specifics  of the calculation we want to perform.
This is done by using ``specwizard.BuildInput``. This class creates the ``Wizard`` python dictionary the is used in all the following steps to calculate spectra.
The class contains a set of functions that create this dictionary, however the easiest way to deal with the input is to feed in a `yml` file and read it with `BuildInput.read_from_yml()`.

In the different examples of usage of ``specwizard`` you will always find the ``yml`` file that produces the output. 
Here we will explain the different elements that are mandatory or optional in the configuration file. 

------------
``file_type``
------------


.. code:: YAML

    file_type:
        
        sim_type: simulation_name
        snap_type: file_type



``sim_type``:  This is the type of simulation from where we will get the data to  construct the sightline. The current supported options are (they are not case sensitive)

- ``eagle``
- ``swift``
- ``colibre``
- ``hydrangea``

``snap_type`` this is the type of file.  There are two options ``snapshot`` and ``LOS``. 

- ``snapshot`` This is the simulation file that contains all the particles from the simulation volume ``specwizard`` will construct a sightline by applying a spatial mask to the particles of the simulation. 
- ``LOS`` This file contains already sightlines pre-made from the simulation output. This method is way faster.

---------------------
``snapshot_paramas`` 
---------------------

this consist of `directory` which is the path that host the simulation data and ``file`` is the name of the .hdf5 file. `specwizard` will perform a check for the existence of both the file and the directory, 

.. code:: YAML

    snapshot_params:
        
        directory: /path/to/the/file
        file: simulation_file.hdf5


-------------
``sightline``
-------------

This section contain all the information that is needed to set-up the sightline. In the case of a full using a  full snapshot we have to define which axis will correspond to the axis of the line of sight.
Example:


 .. code:: YAML

        sightline:
            ProjectionAxes: ['simx','simy','simz']
            ProjectionStart: [0.5,0.5,0]
            ProjectionLength: 1
            SightLength: null
            ProjectionExtend:
                extend: False
                extendfactor: 3
            nsight: 0


``ProjectionAxes``: Creates the correspondence between the simulation axis and the axis of the spectra line. For example, for ``ProjectionAxes=['simz','simy','simz']``, the sight line will - in the coordinates of the simulation, be along the :math:`x-axis`, and go through the point with :math:`(z,y)`.

``ProjectionStart``: are the starting positions for each coordinate as a fraction of the box size, e.g ``x-position = 0.5`` corresponds to half of the boxsize. 

``ProjectionLength``: this is an optional parameter that establish the length of the sightline. If not provided the default value will be the length of the box. 

``ProjectionExtend``: This optional parameter if ``extend: True`` will insert the spectra into an array by ``extendfactor`` times bigger than the sightline length. Is useful to observe dramatic  broadening effects.

``nsight``: it will be used if we use the saving routine to label the line of sight e,g ``nsight=4``-> LOS_4

In case of using a LOS file. The only necessary parameter is ``nsight``. Which will provide the number of LOS that will be read from the LOS file. And will obtain the rest of parameters from the file.


-------------
``ionparams``
-------------
Here we set all the parameters related with the elements and ions that are going to be taken into account for the Spectra. If needed the ionization fraction will be calculated through the interpolation of the ionization table given by ``table_type``.
Example:

 .. code:: YAML

        ionparams:
            table_type: specwizard_cloudy
            iondir: /cosma/home/dp004/dc-aram1/workdir7/pyspecwizard_main/spwtables/HM12/
            fname: null
            ions: [['Hydrogen', 'H I'],['Helium','He II'],['Carbon','C III'],['Carbon','C IV'],['Oxygen','O VI'],['Silicon','Si IV']]
            SFR_properties:
                modify_particle: True
                ignore_particle: False
                Temperature [K]: 1.0E+4    # Putting the plus sign between E and the number is important to not be processed as a string. 
            atomfile: atom_info.hdf5


``table_type``:  Name of the ionizations tables to use to calculate the ion fractions. Current options:

    - ```'specwizard_cloudy'```:  Corresponds to the Cloudy tables created from the Cloudy notebook based in `Haardt & Madau 2012 <https://ui.adsabs.harvard.edu/abs/2012ApJ...746..125H/abstract>`_.
    -  ```'ploeckinger'```: Correspond to the tables described in  `Ploeckinger & Schaye 2020 <https://ui.adsabs.harvard.edu/abs/2020MNRAS.497.4857P/abstract>`_ .

``iondir``: Path for the ionization table.

``fname``: Name of ionization table, only needed for the Ploeckinger tables. 

``ions``: List of elements and Ions to calculate the absorption  e,g ``[[Hydrogen, H I],[Oxygen, O VI],...]``

``SFR_properties``:  if ``modify_particle: True`` we will deal with Star Forming Particles in one of two ways:

 - If ``ignore_particle = True``, will set to zero the IonizationFraction of SFR particles. 
 - If ``ignore_particle = False`` we will set the Temperature of SFR by the value indicated in ``'Temperature [K]'``. 


-------------
``ODParams``
-------------
Here we establish optional parameters to turn on or off effects in the calculation of the optical depth:

Example:

 .. code:: YAML

        ODParams:
            VelOffset_kms: 0
            PecVelEffectsOff: False
            ThermalEffectsOff: False
            VoigtOff: False

``VelOffset_kms``: Generates an velocity offset in kilometers per second that is added to the peculiar velocities of the line.  

``PecVelEffectsOff``: If ``True`` there will be no velocity displacements in the line caused by peculiar velocities. 

``ThermalEffectsOff``: If ``True`` the thermal broadening of the line will be off. 

``VoigtOff``: If ``True`` there will be no damping wings added into the H I line. 

---------------
``LongSpectra``
---------------

This field sets all the necessary parameters to construct a longspectra. Thus is only necessary if you want to compute a LongSpectra. 


 .. code:: YAML

        LongSpectra:
            lambda_min: 300.0
            lambda_max: 8000.0
            dlambda: 0.5
            z_qsr: 3.2
            delta_z: 0.1
            all_contaminants: True
            file_dir: /cosma7/data/Eagle/ScienceRuns/Planck1/L0100N1504/PE/REFERENCE/data/los/ #./Analisis/long_los/ #s

``lambda_min``: Minimal wavelength in Ångströms that will be detected by the "spectrometer".

``lambda_max``: Maximal wavelength in Ångströms that will be detected by the "spectrometer"

``dlambda``: Spectrometer pixel size in Ångströms.

``z_qsr``: Redshift location of the quasar that is used as "backlight" for the LOS. 

``delta_z``: This is the redshift tolerance. Meaning that it will look for all the files that are within current z + delta_z. For more details look at the entry on longspectra. 

``all_contaminants``: If ``True`` will add all the ions available that contribute to the wavelength range. 

``file_dir``: Directory that has all the lines of sight files to construct the longspectra. 