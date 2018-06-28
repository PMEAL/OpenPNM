.. _phys_models_api:

================================================================================
Physics Models
================================================================================

.. contents:: Contents of this Page
    :depth: 3

--------------------------------------------------------------------------------
Overview of Submodule
--------------------------------------------------------------------------------

.. automodule:: openpnm.models.physics

--------------------------------------------------------------------------------
Detailed Model Documentation
--------------------------------------------------------------------------------

................................................................................
Capillary Pressure
................................................................................

.. autofunction:: openpnm.models.physics.capillary_pressure.washburn
.. autofunction:: openpnm.models.physics.capillary_pressure.purcell

----

................................................................................
Diffusive Conductance
................................................................................

.. autofunction:: openpnm.models.physics.diffusive_conductance.ordinary_diffusion

----

................................................................................
Electrical Conductance
................................................................................

.. autofunction:: openpnm.models.physics.electrical_conductance.series_resistors

----

................................................................................
Hydraulic Conductance
................................................................................

.. autofunction:: openpnm.models.physics.hydraulic_conductance.hagen_poiseuille

----

................................................................................
Thermal Conductance
................................................................................

.. autofunction:: openpnm.models.physics.thermal_conductance.series_resistors

----

................................................................................
Multiphase
................................................................................

.. autofunction:: openpnm.models.physics.multiphase.conduit_conductance
.. autofunction:: openpnm.models.physics.multiphase.late_pore_filling
.. autofunction:: openpnm.models.physics.multiphase.late_throat_filling
