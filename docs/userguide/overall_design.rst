.. _overall_design:

###############################################################################
Overall Design
###############################################################################

OpenPNM separates different types of data between 5 object types: **Network**, **Geometry**, **Phase**, **Physics**, and **Algorithms**.  Each of these are described in more detail below, but their names hopefully indicate what sort of data or roles are assigned to each.

The main motivation for this division of data between objects is encompassed by the following table (or grid) and explained below:

+-------------+-----------+-----------+-----------+
| **Network** |  Phase 1  |  Phase 2  |  Phase 3  |
+-------------+-----------+-----------+-----------+
| Geometry 1  | Physics 1 | Physics 2 | Physics 3 |
+-------------+-----------+-----------+-----------+
| Geometry 2  | Physics 4 | Physics 5 | Physics 6 |
+-------------+-----------+-----------+-----------+

This grid represents a single Project.  Each Project has one Network, which has *Np* pores and *Nt* throats.  The Network's main role is to house the pore coordinates and throat connection data.  Because there is only one Network, it occupies the special corner location in the above grid.

A Project can have many Phases, and since each Phase has a different value for a given property (e.g. density or viscosity) a unique object is required for each one.  Each Phase represents a new column in the grid, where each column has unique values of thermo-physical properties.  Phases can exist everywhere, anywhere, or nowhere in a given domain, and can redistribute during a simulation.  As such, Phase properties are calculated everywhere, so they are associated with all pores and throats in the domain.

In some cases, the domain may have multiple distinct regions, such as a two-layered electrode, or multi-modal pore size distributions such as a hierarchical rock.  Since Geometry objects are responsible for calculating the pore and throat sizes, it is necessary to have multiple objects for these cases (e.g. different parameters for the distribution functions).  Each Geometry represents a new row in the grid, where each row has unique values of geometrical properties.  Each row **also** represents a subset of the total pores and throats in the domain, since each pore and throat can only be assigned to one Geometry. Thus Geometry objects have their own values of *Np* and *Nt*, corresponding to the subset of pores and throats they are in charge of.

Finally, Physics objects exist at the intersection of a row and a column.  This represents the fact that a Physics object calculates values that require size information *and* thermo-physical properties.  For example, the Hagan-Poiseuille model for hydraulic conductance requires throat diameter and length, as well as viscosity.  Each Physics object is associated with a specific Phase, whom is retrieves thermo-physical property data from, and a specific Geometry from whom it retries geometrical information.  Physics objects, because they are associated one-to-one with a Geometry, also apply to a subset of pores and throats, hence have their own values of *Np* and *Nt*.

With this grid analogy in mind, we can now dive into an explanation of each object and it's particular abilities.

================================================================================
Object Inheritance Structure
================================================================================

OpenPNM consists of 5 main object types: Network, Phases, Geometries, Physics, and Algorithms.  The inheritance structure of each of these objects is shown in the diagram below.  Each of these objects is a subclass of the :ref:`base_api` class, described in more detail in the next section.  Some of these object also have the ability to store pore-scale models added via the :ref:`modelsmixin_api` mixin class.  Finally, some objects are applied only to subdomains rather than the entire domains, so these inherit from the :ref:`subdomain_api` class, which is itself a subclass of Base.

.. image:: /../docs/static/images/Overall_Inheritance_Diagram.png
    :width: 800px
    :align: center

================================================================================
The Base Class
================================================================================

All the objects in OpenPNM are subclasses of a single :ref:`base_api` class, which is itself a subclass of the `Python Dictionary <http://www.pythonforbeginners.com/dictionary/how-to-use-dictionaries-in-python>`_ (``dict``).  So before explaining each of the specific OpenPNM subclasses, the Base class should be covered.

.. autoclass:: openpnm.core.Base
    :noindex:

================================================================================
Networks
================================================================================

The :ref:`generic_network_api` class add more methods to the Base class than any other type in OpenPNM.  These added methods are all related to the querying of topological information such as finding neighboring throats, or nearby pores. The table below gives a high level overview of these methods.  For a deeper discussion of the topological data format used by OpenPNM (and thus how these queries are performed) refer to :ref:`topology`.

.. autoclass:: openpnm.network.GenericNetwork
    :noindex:

================================================================================
Phases and the ModelsMixin
================================================================================

The :ref:`generic_phase_api` class is very simple subclass of `The Base Class`_.  The subclass itself adds *no* additional methods beyond those of Base, *but* it uses multiple inheritance, so inherits 3 methods from :ref:`modelsmixin_api`, and an added attribute called ``models`` which is a :ref:`modelsdict_api` object that stores the models and their respective parameters.

.. autoclass:: openpnm.phases.GenericPhase
    :noindex:

================================================================================
Subdomains: Geometry and Physics
================================================================================

Geometry and Physics objects are the only two object types can be assigned to a subset of the full domain. This ability is included in the :ref:`subdomain_api` class which is a child of the normal Base class.  The only functionality added to ``Subdomain`` is the ability to set and remove which locations (pores and throats) the object is assigned to.

================================================================================
Algorithms
================================================================================

Like the other classes discussed above, the :ref:`generic_algorithm_api` class inherits from Base, but because every algorithm is a bit different and they tend to be more complicated, the details won't be discussed here.  













.
