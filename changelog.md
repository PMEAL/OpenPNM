# Changes in V2

* All module names, including the package itself, have been changed to lower-case in accord with pep8 rules

* A Project class has been added, that acts as a container for a given pore network simulation, such that the Workspace is populated with a collection of Projects, each of which contains a Network object and all other associated objects (i.e. geometries, phases, etc.)

* pore-scale models can now be assigned to any object, so a model that calculates pore diameters can be assigned to a network if desired

* pore-scale models are now stored in a top level module, making then more prominent and easier to find

* A settings attribute has been added to all objects to contain various flags and values which used to be stored as hidden attributes

* It's now possible to initialize any object independent of supplying existing objects as arguments.

* Added support for two new data export formats: hdf5 and xdmf.  These are designed for storing massive data sets with maximum speed, and the xdmf format can be open in Paraview for vastly sped-up visualization

* Topological manipulation tools have been moved from a sub-module of network to a top level module (topotools) to make it easier to access and more prominent

*
