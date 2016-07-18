# CHANGE LOG

## Version 1.5 (June 1, 2016)

### Major Enhancements and Changes

1. Added two new import classes for importing from [iMorph](http://imorph.fr) and [3DMA-Rock](http://www.ams.sunysb.edu/~lindquis/3dma/3dma_rock/3dma_rock.html).

2. Now uses [dill](https://pypi.python.org/pypi/dill) instead of Python's standard [pickle](https://docs.python.org/3.5/library/pickle.html) library.  This means that custom pore scale models and classes can be saved to "pnm" files.

3. Changed the name of the "Controller" object to "Workspace", which more literally describes it's role.  "Controller" can still be used for backward compatibility.

4. All Network and topology manipulation tools are now found under ```Network.tools``` although they can still be accessed via ```Utilities.topology``` for backwards compatibility.

5. Completely reworked documentation.  The 'user guide' now consists of 3 tutorials of varying levels, with the aim that users will know how to use OpenPNM pretty well after completing all three.  There is also a separate 'reference' section in the user guide that explains the inner workings of the code in more detail, but users don't have to sift through these to get started anymore.

6. All documentation is now hosted on [ReadTheDocs](https://readthedocs.org/projects/openpnm/), which rebuilds the documentation every time a new commit is pushed to Github.  The docs will now always be up-to-date!

7. A new parallel repository has been created to house all Examples at https://github.com/PMEAL/OpenPNM-Examples.  All the code in this repository is tested against the latest version of OpenPNM on PyPI, so if there are any broken examples we'll know about it.  This will remedy the frustration felt by many users when trying to learn OpenPNM by examples.  This is also why the Tutorials were created as the main User Guide (see 5 above).

### New Features

1. Added ability to have spatially varying pore seed values, which is useful for creating porosity distributions or undulations.

2. Added new arguments to ```num_neighbors```.  It can now apply set logic (i.e. 'union' and 'intersection') when counting neighbors, and it can now count neighboring throats as well as pores with the ```element``` keyword (which is 'pores' by default to maintain backwards compatibility).

3. Object handles are now stored in dicts rather than lists, so you can access them by name (geom.phases['air']) and also iterate on the dict (pn.phases.values()). This was implemented in a backwards compatible way so geom.phases() and geom.phases('air') still work as they used to...of course using the dict syntax is encourage henceforth.

### Minor Improvements and Bug Fixes

1. Fixed code coverage reporting on [Codecov](https://codecov.io/gh/PMEAL/OpenPNM)

2. Changed many print statements to logger messages

3. Several private methods were removed that were never called by the code

### Acknowledgements

The OpenPNM Developers would like to thank the following people for contributing to this release:

Matthew Stadelman (@stadelmanma) for crafting the iMorph IO class

Masa Prodanovic for help with the 3DMA-Rock import class

## Version 1.4 (February 22nd, 2016)

### Major Enhancements and Additions

1. Added a new Drainage algorithm with much more functionality than *OrdinaryPercolation*.  This algorithm was completely overhauled and designed specifically for simulating drainage experiments.

2. Major upgrade to the Network Importing and Export capabilities, including the ability to export data to a Pandas *DataFrame*, and import from NetworkX.

### New Features

1. ``merge_pores`` was added to *Utilities.topology*

2. ``add_boundary_pores`` and ``make_periodic_connections`` were added the Cubic class to provide more power over the process of adding boundary pores compared to the ``add_boundaries`` method.

3.  ``find_path`` was added to the *Utilities.misc* for finding paths between pairs of pores.

4.  Added new ``pore_diameter`` models named according to the specific distribution they produce, which is much easier to use than the previous general functions.

5.  Added an Empty network class which accepts Nt and Np arguments, so can be used to manually build a network.

6.  Enhanced the ``clear`` method on the Core class to provide more control when clearing data from an object.

### Minor Improvements

1. Improvements to the performance of some method in *GenericNetwork*

2. Various bug fixes in the ``subdivide`` method and *GenericLinearTransport*

3. Added several private parse methods to the *Core* class that help validate inputs to the various methods and provide more helpful error message.

4.  Altered the *GenericPhysics* class so that its associated *Phase* can be changed.

## Version 1.3 (July 9th, 2015)

### Major Enhancements and Additions

1. The main upgrades have been in the *Delaunay* network generation class, and the related *Voronoi* geometry class and related models. These changes center around the use of image analysis to determine pore and throat sizes, and creating a voxelized representation of the solid structure for visualization. A detailed example has been created on the use of the Delaunay-Voronoi class and will be posted to http://openpnm.org.

2. Additionally, this version now has much more extensive test coverage exceeding 80%.

## Version 1.2 (July 9th, 2015)

### Major Enhancements and Additions

1. Major Reorganization of the *Network* topology manipulations methods. There is now a *topology* class in the *Utilities* submodule, which houses all manipulation methods such as extend, trim, stitch and so on. Helper/wrapper methods are present on Network objects to maintain backwards compatibility.

2. Vastly improved test coverage in the form of adding a unit test framework, TravisCI integration, Coveralls coverage checking, PEP8 compatibility

### New Features

1. * Added ``find_nearest_pores`` to GenericNetwork for finding spatially nearby pores regardless of whether they are topologically connected.

2. Added a ``subdivide`` method to the suite of topology manipulation tools, which allows single pores to be divided into many pores to produce hierarchical networks.

3. New Network methods: ``find_nearby_pores`` and ``find_clusters2``

## Version 1.1 (April 9th, 2015)

* Improved the interaction with models on each object.
* Introduced the *Controller* class, which provides high level oversight to all simulations.  The *Controller* allows saving and loading of simulations and entire workspaces.
* Added the ability apply source and sink terms to the transport solvers.  The ``set_source_term`` method was added to *GenericLinearTransport* class and a variety of generic source term models have been added to the available Physics models.
