API reference
=============

Pipeline
--------

The functions and classes most users need, in the order they are used.

.. currentmodule:: specwizard

.. autosummary::
   :toctree: generated
   :nosignatures:

   Build_Input
   GenerateShortSpectra
   LongSpectra
   OpticalDepth_IO
   Analyse_Opticaldepth

Lower-level building blocks
---------------------------

.. autosummary::
   :toctree: generated
   :nosignatures:

   ReadData
   SightLineProjection
   ComputeOpticaldepth
   IonTables
   Elements
   Atomfile
   Run_cloudy
   read_obs_data
   ReadPhys

Internals
---------

Classes used inside the pipeline; useful if you want to compute line profiles
by hand or extend SPECWIZARD to a new simulation code.

.. autosummary::
   :toctree: generated
   :nosignatures:

   SpecWizard_Lines.Lines
   h1damping.HIDamping
   SpecWizard_IonizationBalance.IonizationBalance
   SpecWizard_SplineInterpolation.Bspline
   SpecWizard_SplineInterpolation.TGauss
   SpecWizard_SplineInterpolation.ColumnTable
   reading_simulations.InputFunctions
   reading_simulations.ReadEagle
   reading_simulations.ReadSwift
   reading_simulations.ReadHydrangea
   SimulationInputKeys.get_simkeys
