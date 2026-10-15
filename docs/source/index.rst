SPECWIZARD
==========

.. include:: ../../README.rst
   :start-after: INTRO_FLAG
   :end-before: INTRO_FLAG_END

.. grid:: 1 2 2 2
   :gutter: 3

   .. grid-item-card:: Getting started
      :link: getting_started/installation
      :link-type: doc

      Install SPECWIZARD, get the atomic data and ionization tables, and run
      your first spectrum.

   .. grid-item-card:: User guide
      :link: user_guide/configuration
      :link-type: doc

      Every section of the ``Wizard.yml`` configuration file, and the layout of
      the HDF5 files SPECWIZARD writes.

   .. grid-item-card:: Tutorials
      :link: tutorials/index
      :link-type: doc

      Worked Jupyter notebooks: short spectra, building an atom file, Lyman-limit
      systems.

   .. grid-item-card:: API reference
      :link: api/index
      :link-type: doc

      Auto-generated documentation of every public class and function.

.. toctree::
   :hidden:
   :caption: Getting started

   getting_started/installation
   getting_started/quickstart

.. toctree::
   :hidden:
   :caption: User guide

   user_guide/configuration
   user_guide/output_format

.. toctree::
   :hidden:
   :caption: Tutorials

   tutorials/index

.. toctree::
   :hidden:
   :caption: Reference

   api/index
