Overview
========
.. INTRO_FLAG

SPECWIZARD is a Python package to compute and analyze mock quasar absorption spectra from cosmological simulations.

This repository provides:

- Tools to build and validate simulation input dictionaries ("Wizard" dictionaries).
- Routines to generate short spectra for single sightlines.
- Routines to assemble full long spectra across a redshift path.
- HDF5 I/O helpers to save and read spectra in a consistent format.


.. INTRO_FLAG_END

Citing SPECWIZARD
-----------------

If you use SPECWIZARD in scientific work, please cite the relevant code/paper for your release.


Installation Notes
==================

.. INSTALL_FLAG

SPECWIZARD is written in ``python3`` and stable versions are available on PyPI_.
The easiest installation method is:


.. code-block::

    pip install specwizard 

The code is constantly updated. Feedback and bug reports are welcome through email_ or GitHub.

.. _PyPI: https://pypi.org/project/specwizard/
.. _email: mailto:aramburo@lorentz.leidenuniv.nl 
.. INSTALL_FLAG_END


Quick Start
===========

Most users follow this flow:

1. Create a YAML configuration file.
2. Load it with ``Build_Input.read_from_yml``.
3. Generate spectra.
4. Save to HDF5.
5. Read the file back for analysis.


Minimal Configuration (Wizard YAML)
===================================

Save a file such as ``Wizard.yml``:

.. code-block:: yaml

        file_type:
            sim_type: swift
            snap_type: snapshot

        snapshot_params:
            directory: /path/to/snapshots
            file: snapshot_012.hdf5

        sightline:
            ProjectionAxes: [simx, simy, simz]
            ProjectionStart: [0.5, 0.5, 0.0]
            ProjectionLength: 1
            SightLength: null
            ProjectionExtend:
                extend: false
                extendfactor: 3
            nsight: 0

        ionparams:
            table_type: specwizard_cloudy
            iondir: /path/to/ion_tables
            fname: null
            ions:
                - [Hydrogen, H I]
                - [Carbon, C IV]
            SFR_properties:
                modify_particle: true
                ignore_particle: false
                Temperature [K]: 1.0E+4
            atomfile: /path/to/atomic_info.hdf5

        ODParams:
            VelOffset_kms: 0
            PecVelEffectsOff: false
            ThermalEffectsOff: false
            VoigtOff: false

        extraparams:
            periodic: true
            pixkms: 1
            ReadIonFrac:
                ReadIonFrac: false
                ReadHydrogen: true
                HI: NeutralHydrogenAbundance
                ReadHelium: false
                He: ""
                fname_urchin: ""

        Output:
            directory: ./outputs/
            fname: shortspec.hdf5

Notes:

- ``Build_Input.read_from_yml`` reads ``file_type``, ``snapshot_params``, ``sightline``, ``ionparams``, ``ODParams`` and optional ``LongSpectra``, ``extraparams``, and ``Output`` sections.
- ``Output`` is required if you plan to use ``OpticalDepth_IO`` to save data.


Generate And Save A Short Spectrum
==================================

.. code-block:: python

        import specwizard as spw

        # 1) Build Wizard dictionary from YAML
        builder = spw.Build_Input()
        wizard = builder.read_from_yml("Wizard.yml")

        # 2) Generate short-spectrum products
        opticaldepth, projected_los, particles = spw.GenerateShortSpectra(wizard)

        # 3) Save to HDF5 using the package writer
        io = spw.OpticalDepth_IO(wizard=wizard, create=True)
        payload = {
                "nsight": wizard["sightline"]["nsight"],
                "Projection": projected_los,
                "OpticaldepthWeighted": opticaldepth,
        }
        io.write_shortspectra_to_file(payload)


Read A Saved Short Spectrum
===========================

.. code-block:: python

        import specwizard as spw

        builder = spw.Build_Input()
        wizard = builder.read_from_yml("Wizard.yml")

        io = spw.OpticalDepth_IO(wizard=wizard, create=False)
        short_data = io.read_shortspectra_from_file()

        # Examples of access:
        header = short_data["Header"]
        one_los = next(iter(short_data["Data"]))
        los_data = short_data["Data"][one_los]


Generate, Save, And Read A Long Spectrum
========================================

To generate a long spectrum, include a ``LongSpectra`` block in your YAML:

.. code-block:: yaml

        LongSpectra:
            lambda_min: 945.0
            lambda_max: 8000.0
            dlambda: 0.5
            z_qsr: 3.0
            delta_z: 0.01
            all_contaminants: false
            file_dir: /path/to/los/files/

Then run:

.. code-block:: python

        import specwizard as spw

        builder = spw.Build_Input()
        wizard = builder.read_from_yml("Wizard.yml")

        long_builder = spw.LongSpectra(wizard)
        coven, redshifts = long_builder.create_coven()
        long_spectra = long_builder.do_long_spectra(coven)

        # Optional post-processing
        # long_spectra = long_builder.add_contaminants(long_spectra)
        # long_spectra = long_builder.add_HI_damping_wings(long_spectra, n=2)

        # Save
        io = spw.OpticalDepth_IO(wizard=wizard, create=True)
        io.write_fullspectrum_to_file(long_spectra)

        # Read
        io_read = spw.OpticalDepth_IO(wizard=wizard, create=False)
        long_data = io_read.read_fullspectrum_from_file()


Data Fields: What You Can Read
================================

This project writes structured HDF5 groups for both short and long spectra. Below are the common fields and where to find them.

- Short spectra (per LOS):
    - Top-level groups: ``LOS_<n>/`` and ``Header/``.
    - Element-weighted properties: ``LOS_<n>/<Element>/Element-weighted/<field>`` (e.g. ``Velocities``, ``Densities``, ``Temperatures``).
    - Ion-weighted properties: ``LOS_<n>/<Element>/<Ion>/Ion-weighted/<field>``.
    - Optical-depth-weighted properties for each ion: ``LOS_<n>/<Element>/<Ion>/Ion optical depth-weighted/<field>`` (includes ``Optical depths``, ``Velocities``, ``Densities``, ``Temperatures``, ``Metallicities``, ``HydrogenDensities`` where present).
    - If the simulation includes non-equilibrium / tracked ionic abundances, additional groups are written under ``SimIon-weighted`` and ``SimIon optical depth-weighted`` with the same fields as the tabulated ions.

- Long spectra (full-spectrum accumulation):
    - Global grids: ``FullSpectrum/Velocities`` and ``FullSpectrum/Wavelengths``.
    - Per-ion data: ``FullSpectrum/<Element>/<Ion>/...`` with datasets for ``lambda0``, ``f-value`` and fields such as ``Optical depths``, ``Velocities``, ``Densities``, ``Temperatures``, ``Metallicities``, ``HydrogenDensities`` (each saved together with metadata/attributes when available).
    - The long-spectrum writer prefers simulation-tracked ion optical-depth-weighted fields when present: if ``SimIons`` optical-depth-weighted quantities exist they will be used (and saved) in preference to tabulated ion optical-depth-weighted values.

- Access helpers:
    - Use ``OpticalDepth_IO.ReadVariable(path)`` to read a specific dataset and its attributes.
    - ``OpticalDepth_IO.ReadHeader()`` returns header attributes and any stored global datasets.

Output And File Layout
======================

- Output path is controlled by ``wizard['Output']['directory']`` and ``wizard['Output']['fname']``.
- Short-spectrum files contain LOS-based groups such as ``LOS_0/...`` plus ``Header/...``.
- Full long-spectrum files are stored under ``FullSpectrum/...``.
- ``OpticalDepth_IO.ReadVariable(path)`` can be used to load a specific dataset with metadata.

