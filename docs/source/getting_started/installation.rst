Installation
============

.. include:: ../../../README.rst
   :start-after: INSTALL_FLAG
   :end-before: INSTALL_FLAG_END

To install the development version from GitHub instead:

.. code-block:: bash

   git clone https://github.com/specwizard/specwizard.git
   cd specwizard
   pip install -e .

Required data files
-------------------

Besides the simulation output itself, SPECWIZARD needs two auxiliary files.
Their paths are set in the ``ionparams`` section of the configuration file
(see :doc:`../user_guide/configuration`).

Atomic data file (``atomfile``)
   An HDF5 file with the rest wavelength :math:`\lambda_0` and oscillator
   strength :math:`f` of every transition. You can download the pre-built file
   or generate your own from NIST with :class:`specwizard.Atomfile`; see the
   :doc:`../tutorials/atomic_data` tutorial.

   .. TODO: add the download link for the pre-built atomic_info.hdf5

Ionization tables (``iondir`` / ``fname``)
   Tables of ion fractions as a function of density, temperature and redshift.
   Supported types are ``specwizard_cloudy`` (Haardt & Madau 2012),
   ``ploeckinger`` (Ploeckinger & Schaye 2020) and ``chimes``
   (Ploeckinger & Richings 2025). The ``specwizard_cloudy`` tables can be
   generated with :class:`specwizard.Run_cloudy`.

   .. TODO: add download links for the ionization tables

Building this documentation locally
-----------------------------------

.. code-block:: bash

   pip install -r docs/requirements.txt
   make -C docs html          # output in docs/_build/html/index.html
