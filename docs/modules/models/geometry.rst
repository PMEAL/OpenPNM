.. _geom_models_api:

================================================================================
Geometry Models
================================================================================

.. contents:: Contents of this Page
    :depth: 3

--------------------------------------------------------------------------------
Overview of Submodule
--------------------------------------------------------------------------------

.. automodule:: openpnm.models.geometry

--------------------------------------------------------------------------------
Detailed Model Documentation
--------------------------------------------------------------------------------

................................................................................
Pore Size
................................................................................

.. autofunction:: openpnm.models.geometry.pore_size.weibull
.. autofunction:: openpnm.models.geometry.pore_size.normal
.. autofunction:: openpnm.models.geometry.pore_size.generic
.. autofunction:: openpnm.models.geometry.pore_size.largest_sphere
.. autofunction:: openpnm.models.geometry.pore_size.equivalent_diameter

----

................................................................................
Pore Volume
................................................................................

.. autofunction:: openpnm.models.geometry.pore_volume.sphere
.. autofunction:: openpnm.models.geometry.pore_volume.cube

----

................................................................................
Pore Seed
................................................................................

.. autofunction:: openpnm.models.geometry.pore_seed.random
.. autofunction:: openpnm.models.geometry.pore_seed.spatially_correlated

----

................................................................................
Pore Surface Area
................................................................................

.. autofunction:: openpnm.models.geometry.pore_surface_area.sphere
.. autofunction:: openpnm.models.geometry.pore_surface_area.cube

----

................................................................................
Throat Area
................................................................................

.. autofunction:: openpnm.models.geometry.throat_area.cylinder
.. autofunction:: openpnm.models.geometry.throat_area.cuboid

----

................................................................................
Throat Size
................................................................................

.. autofunction:: openpnm.models.geometry.throat_size.weibull
.. autofunction:: openpnm.models.geometry.throat_size.normal
.. autofunction:: openpnm.models.geometry.throat_size.generic
.. autofunction:: openpnm.models.geometry.throat_size.from_neighbor_pores
.. autofunction:: openpnm.models.geometry.throat_size.equivalent_diameter

----

................................................................................
Throat Length
................................................................................

.. autofunction:: openpnm.models.geometry.throat_length.piecewise

----

................................................................................
Throat Perimeter
................................................................................

.. autofunction:: openpnm.models.geometry.throat_perimeter.cylinder
.. autofunction:: openpnm.models.geometry.throat_perimeter.cuboid

----

................................................................................
Throat Surface Area
................................................................................

.. autofunction:: openpnm.models.geometry.throat_surface_area.cylinder
.. autofunction:: openpnm.models.geometry.throat_surface_area.cuboid
.. autofunction:: openpnm.models.geometry.throat_surface_area.extrusion

----

................................................................................
Throat Volume
................................................................................

.. autofunction:: openpnm.models.geometry.throat_volume.cylinder
.. autofunction:: openpnm.models.geometry.throat_volume.cuboid
.. autofunction:: openpnm.models.geometry.throat_volume.extrusion

----

................................................................................
Throat Shape Factor
................................................................................

.. autofunction:: openpnm.models.geometry.throat_shape_factor.compactness
.. autofunction:: openpnm.models.geometry.throat_shape_factor.mason_morrow
.. autofunction:: openpnm.models.geometry.throat_shape_factor.jenkins_rao
