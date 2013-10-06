===============================================================================
Network Geometry
===============================================================================
The Geometry module contains methods for:
1. Generating networks of various topologies (such as Cubic and Delaunay)
2. Importing networks generated using other code (such as extraction from 3D tomographic images)
3. Altering existing the geometry networks (such as merging networks)

-------------------------------------------------------------------------------
Generating Geometry
-------------------------------------------------------------------------------

*******************************************************************************
Generic Geometry
*******************************************************************************
The GenericGeometry class takes advantage of the object-oriented nature of Python by using `inheritance <http://docs.python.org/2/tutorial/classes.html>`_.  There are three main components to inheritance as used here.  

Firstly, the GenericGeometry class contains methods that will likely be used by all methods regardless of topology.  Every network inherits these methods and *can* use them if desired.  For instance, generate_pore_seeds can be used "as is" to assign a random number to each pore for subsequent use in pore size distribution calculations (i.e. generate_pore_diameters).  Alternatively, if the user wishes to use a more complex method to generate random seeds (for instance with spatial correlation), then they are free to over-write or sub-class this method in their specific generator.  The procedure for accomplishing this is outlined Writing Custom Geometry section below.  

The second component to inheritance as applied in OpenPNM is that some methods *must* be sub-classed.  Specifically, generate_pores() and generate_throats() are not implimented in the GenericGeometry class.  One of the main results of generate_pores is to dictate the spatial coordinates of the pore centers; generate_throats dictates which pores are connected to which neighbors.  A cubic network places pores in space and defines connections in a completely different way than a random network.  Thus, there is no generic or default way to generate different networks.  

And finally, all sub-classed methods have the same name and black-box functionality as the generic methods.  The enables a truly generic generation scheme (defined by the generate method in GenericGeometry) that always calls the same methods but acheives different results depending on which methods have been overwritten.  It should be pointed out however, that even the generate method *can* be sub-classed if desired, so there is no need to adhere to the steps or order defined there.  

The generate() method provided in GenericGeometry contains the following methods in the order given below.  The first 3 of these methods are abstract or empty methods that *must* be subclassed by the specific Geometry class.  The remainder are actually implimented in the GenericGeometry class and perform the most basic or common version of their specific function.  Any of these can also be subclassed in each specific Geometry class. Moreover, if the methods or the order given below are unsuitable the the generate() method itself can be subclassed by a specific Geometry class.  

.. automethod:: OpenPNM.Geometry.GenericGeometry.generate_pores()

.. automethod:: OpenPNM.Geometry.GenericGeometry.generate_throats()

.. automethod:: OpenPNM.Geometry.GenericGeometry.add_boundaries()

.. automethod:: OpenPNM.Geometry.GenericGeometry.generate_pore_seeds()

.. automethod:: OpenPNM.Geometry.GenericGeometry.generate_throat_seeds()

.. automethod:: OpenPNM.Geometry.GenericGeometry.generate_pore_diameters()

.. automethod:: OpenPNM.Geometry.GenericGeometry.generate_throat_diameters()

.. automethod:: OpenPNM.Geometry.GenericGeometry.calc_pore_volumes()

.. automethod:: OpenPNM.Geometry.GenericGeometry.calc_throat_lengths()

.. automethod:: OpenPNM.Geometry.GenericGeometry.calc_throat_volumes()

*******************************************************************************
Cubic
*******************************************************************************
The most common and basic type of pore network is based on cubic geometry, with cubic lattice-type connectivity between pores.  

The BasicCubic geometry corresponds to simplest `Bravais Lattice <http://en.wikipedia.org/wiki/Bravais_lattice>`_ type, the Primitive Centered Cubic, pcc.


*******************************************************************************
Template
*******************************************************************************
This is a varient of the Cubic network that allows for arbitrarily complex shapes such as spheres and cylinders, but still defines connections between pores based on lattice-type connectivity.  

There are two main motivations for including this generator.  Firstly, it is the most straightforward way to generate unusual custom geometry of any shape.  Modeling the coking of catalyst particles of spherical or cylindrical shape can be accomplished with equal ease.  Secondly, some users will be more comfortable dealing with numerical matrices outside of OpenPNM and this generator allows them to store network data in a more human-friendly manner (i.e. in a series of matrices the same shape as the network).  For instance, it is possible to generate cubic networks this way if an image of a cube is provided.  

The ArbitraryCubic geometry generator accepts a binary 3D (or 2D) image with some pattern of 1's to define the network shape.  Generating a spherical network using this generator is accomplished using the ndimage package in Scipy as follows:

.. code-block:: python
     
   import scipy as sp
   import scipy.ndimage as spim
   img = sp.ones((20,20,20),dtype=boolean)
   img[10,10,10] = 0
   img = spim.distance_transform_edt(img)
   img = img<=5
   pn = OpenPNM.Geometry.Custom(image_shape=img).generate()
   
This will generate a spherical network with cubic-lattice connectivity.  All pore and throat properties will be generated from the methods inherited from GenericGeometry.  It is possible to specify custom properties to overwrite those produced by the generic methods.  For instance, if pore sizes are larger near the surface than near the core of the sphere this can be calculated externally, stored in an image of the desired shape, and then imported into the network as follows:

.. code-block:: python
     
   import scipy as sp
   import scipy.ndimage as spim
   img = sp.ones((20,20,20),dtype=boolean)
   img[10,10,10] = 0
   img = spim.distance_transform_edt(img)
   img = img<=5
   pn = OpenPNM.Geometry.Custom(image_shape=img).generate()


*******************************************************************************
Delaunay
*******************************************************************************
This a basic type of random network generated by placing the specified number of basepoints randomly in the domain, and then determine which pores are neighbors using a Delaunay tessellation.  




-------------------------------------------------------------------------------
Importing Networks
-------------------------------------------------------------------------------





-------------------------------------------------------------------------------
Writing Custom Generators
-------------------------------------------------------------------------------

*******************************************************************************
Sub-classing Methods in GenericGeometry
*******************************************************************************


*******************************************************************************
Adding New Methods using generate_misc()
*******************************************************************************




-------------------------------------------------------------------------------
Manipulating Geometry
-------------------------------------------------------------------------------


*******************************************************************************
Extract Sub-Network
*******************************************************************************




*******************************************************************************
Translate and Rotate Network
*******************************************************************************




*******************************************************************************
Stitch Networks
*******************************************************************************





