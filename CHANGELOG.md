# CHANGE LOG

## Version 1.4 (February 22nd, 2016)

### Major Enchancements and Additions

1. Added a new Drainage algorithm with much more functionality than OrdinaryPercolation.  This algorithm was completely overhauled and designed specfically for simulating drainage experiments.

2. Major upgrade to the Network Importing and Export capabilities, including the ability to export data to a Pandas DataFrame, and import from NetworkX.

### New Features

1. ``merge_pores`` was added to Utilities.topology

2. ``add_boundary_pores`` and ``make_periodic_connections` were added the Cubic class to provide more power over the process of adding boundary pores compared to the ``add_boundaries`` method

3.  ``find_path`` was added to the Utilities.misc for finding paths between pairs of pores.

4.  Added an Empty network class which accepts Nt and Np arguments, so can be used to manaully build a network.

5.  Enhanced the ``clear`` method on the Core class to provide more control when clearing data from an object.

### Minor Improvements

1. Improvements to the performance of some method in GenericNetwork

2. Various bug fixes in the subdivide method and the linear solver

3. Added several private 'parse' methods to the Core class that help validate inputs to the various methods and provide more helpful error message.

4.  Altered the GenericPhysics class so that its associated Phase can be changed.

## Version 1.3 (July 9th, 2015)

### Major Enchancements and Additions

1. The main upgrades have been in the Delaunay Network generation class, and the related Voronoi Geometry class and related models. These changes center around the use of image analysis to determine pore and throat sizes, and creating a voxelized representation of the solid structure for visualization. A detailed example has been created on the use of the Delaunay-Voronoi class and will be posted to http://openpnm.org.

2. Additionally, this version now has much more extensive test coverage exceeding 80%.

## Version 1.2 (July 9th, 2015)

### Major Enchancements and Additions

1. Major Reorganization of the Network topology manipulations methods. There is now a topology class in the Utilities submodule, which houses all manipulation methods such as extend, trim, stitch and so on. Helper/wrapper methods are present on Network objects to maintain backwards compatibility.

2. Vastly improved test coverage in the form of adding a unit test framework, TravisCI integration, Coveralls coverage checking, PEP8 compatibility

### New Features

1. * Added ``find_nearest_pores`` to GenericNetwork for finding spatially nearby pores regardless of whether they are topologically connected.

2. Added a subdivide method to the suite of topology manipulation tools, which allows single pores to be divided into many pores to produce hierarchical networks.

3. New Network methods: find_nearby_pores and find_clusters2

## Version 1.1 (April 9th, 2015)

* Improved the interaction with models on each object.
* Introduced the Controller class, which provides high level oversight to all simulations.  The Controller allows saving and loading of simulations and entire workspaces.
* Added the ability apply source and sink terms to the transport solvers.  The ``set_source_term`` method was added to GenericLinearTransport class and a variety of generic source term have been added to the available Physics models.


