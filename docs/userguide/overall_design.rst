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
The Base Class
================================================================================

All the objects in OpenPNM are subclasses of a single :ref:`base_api` class, which is itself a subclass of the `Python Dictionary <http://www.pythonforbeginners.com/dictionary/how-to-use-dictionaries-in-python>`_ (``dict``).  So before explaining each of the specific OpenPNM subclasses, the Base class should be covered.

.. autoclass:: openpnm.core.Base
    :noindex:

================================================================================
Network
================================================================================

The :ref:`generic_network_api` class add more methods to the Base class than any other type in OpenPNM.  These added methods are all related to the querying of topological information such as finding neighboring throats, or nearby pores. The table below gives a high level overview of these methods.  For a deeper discussion of the topological data format used by OpenPNM (and thus how these queries are performed) refer to :ref:`topology`.

.. autoclass:: openpnm.network.GenericNetwork
    :noindex:





.
