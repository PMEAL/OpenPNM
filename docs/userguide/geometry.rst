*******************************************************************************
Network Geometry
*******************************************************************************
The Geometry module contains methods for:
1. Generating networks of various topologies (such as Cubic and Delaunay)
2. Importing networks generated using other code (such as extraction from 3D tomographic images)
3. Altering existing the geometry networks (such as merging networks)

===============================================================================
Generating Geometry
===============================================================================

-------------------------------------------------------------------------------
Generic Geometry
-------------------------------------------------------------------------------
The GenericGeometry class takes advantage of the object-oriented nature of Python by using `inheritance <http://docs.python.org/2/tutorial/classes.html>`_.  There are three main components to inheritance as used here.  

Firstly, the GenericGeometry class contains methods that will likely be used by all methods regardless of topology.  Every network inherits these methods and *can* use them if desired.  For instance, generate_pore_seeds can be used "as is" to assign a random number to each pore for subsequent use in pore size distribution calculations (i.e. generate_pore_diameters).  Alternatively, if the user wishes to use a more complex method to generate random seeds (for instance with spatial correlation), then they are free to over-write or sub-class this method in their specific geometry class.  The procedure for accomplishing this is outlined in the Writing Custom Geometry section below.  

The second component to inheritance as applied in OpenPNM is that some methods *must* be sub-classed.  Specifically, generate_pores() and generate_throats() are not implimented in GenericGeometry.  One of the main results of generate_pores is to dictate the spatial coordinates of the pore centers; generate_throats dictates which pores are connected to which neighbors.  A cubic network places pores in space and defines connections in a completely different way than a random network.  Thus, there is no generic or default way to generate different networks.  

And finally, all sub-classed methods have the same name and black-box functionality as the generic methods.  The enables a truly generic generation scheme (defined by the generate() method in GenericGeometry) that always calls the same methods but acheives different results depending on which methods have been sub-classed.  It should be pointed out however, that even the generate() method *can* be sub-classed if desired, so there is no need to adhere to the steps or order defined there.  

The generate() method provided in GenericGeometry contains the following methods in the order given below.  The first 4 of these methods are abstract or empty methods that *must* be subclassed by the specific Geometry class.  The remainder are actually implimented in the GenericGeometry class and perform the most basic or common version of their specific function.  Any of these can also be subclassed in each specific Geometry class. Moreover, if the methods or the order given below are unsuitable then the generate() method itself can be subclassed by a specific Geometry class.

.. automethod:: OpenPNM.Geometry.GenericGeometry._generate_setup

.. automethod:: OpenPNM.Geometry.GenericGeometry._generate_pores

.. automethod:: OpenPNM.Geometry.GenericGeometry._generate_throats

.. automethod:: OpenPNM.Geometry.GenericGeometry._add_boundaries

.. automethod:: OpenPNM.Geometry.GenericGeometry._generate_pore_seeds

.. automethod:: OpenPNM.Geometry.GenericGeometry._generate_throat_seeds

.. automethod:: OpenPNM.Geometry.GenericGeometry._generate_pore_diameters

.. automethod:: OpenPNM.Geometry.GenericGeometry._generate_throat_diameters

.. automethod:: OpenPNM.Geometry.GenericGeometry._calc_pore_volumes

.. automethod:: OpenPNM.Geometry.GenericGeometry._calc_throat_lengths

.. automethod:: OpenPNM.Geometry.GenericGeometry._calc_throat_volumes

-------------------------------------------------------------------------------
Cubic
-------------------------------------------------------------------------------
The most common and basic type of pore network is based on cubic geometry, with cubic lattice-type connectivity between pores.  The Cubic geometry corresponds to simplest `Bravais Lattice <http://en.wikipedia.org/wiki/Bravais_lattice>`_ type, the Primitive Centered Cubic, pcc.  Each pore is connected to 6 neighbors (in 3D).

-------------------------------------------------------------------------------
Template
-------------------------------------------------------------------------------
This is a varient of the Cubic network that allows for arbitrarily complex domain shapes such as spheres and cylinders, but still defines connections between pores based on lattice-type connectivity.  

There are two main motivations for including this generator.  Firstly, it is the most straightforward way to generate unusual geometry of any shape.  Modeling the coking of catalyst particles of spherical or cylindrical shape can be accomplished with equal ease.  Secondly, some users will be more comfortable dealing with numerical matrices outside of OpenPNM and this generator allows them to store network data in a more human-friendly manner (i.e. in a series of matrices the same shape as the network).  For instance, it is possible to generate cubic networks this way if an image of a cube is provided.  

The Template geometry generator accepts a 3D or 2D ndarray with some pattern of 1's to define the network shape.  Generating a spherical network using this generator can be accomplished using the ndimage package in Scipy as follows:

.. code-block:: python
     
   import scipy as sp
   import scipy.ndimage as spim
   sphere = sp.ones((21,21,21),dtype=boolean)
   sphere[11,11,11] = 0
   sphere = spim.distance_transform_edt(sphere)
   template = sphere<=5.0
   params = {'template' = template}
   pn = OpenPNM.Geometry.Template().generate(**params)
   
This will generate a spherical network with cubic-lattice connectivity.  All pore and throat properties will be generated from the methods inherited from GenericGeometry.  It is possible to specify certain properties in place of or in addition to those produced by the Generic methods.  For instance, if pore sizes are larger near the surface than near the core of the sphere this can be calculated externally, stored in an ndarray of the desired shape, and then imported into the network as follows:

.. code-block:: python

   radial_position = (sphere/5.0)**(0.2)
   pdia_template = radial_position*params{'template'}
   OpenPNM.Geometry.Template().add_pore_property_from_template(pn,pdia_template,'diameter')


-------------------------------------------------------------------------------
Delaunay
-------------------------------------------------------------------------------
This a type of random network generated by placing the specified number of basepoints randomly in the domain, and then determining which pores are neighbors using a Delaunay tessellation.  



===============================================================================
Importing Networks
===============================================================================
...
..
.

===============================================================================
Customizing Existing Geometries
===============================================================================

-------------------------------------------------------------------------------
Sub-classing Methods in GenericGeometry
-------------------------------------------------------------------------------
The ability to subclass methods from a generic class enables very simple customization.  To illustrate the process of sub-classing, let's say we wish to calculate pore volumes assuming they are cubes rather than spheres (which is the default behavior in GenericGeometry) and let's assuming say we want to apply this to the Cubic geometry, but none of the others.  

We begin by noting that pore volumes are calcuated by the _calc_pore_volumes() method in GenericGeometry.  We also note that this method is called during the GenericGeometry._generate() stage.  We do not wish to change the generic behavior for volume calculation or generation, only the behavior of the Cubic geometry. Accordingly we add a method to the Cubic geometry file called _calc_pore_volumes() where we can define the desired volume calculation equations.  It will look something like this:

.. code-block:: python
   def _calc_pore_volumes(self):
       self._net.pore_properties['volume'] = self._net.pore_properties['diameter']**3

When the program is executed, the version of _calc_pore_volumes() located in Cubic will be run rather than the one in GenericGenerator.  

-------------------------------------------------------------------------------
Adding New Methods
-------------------------------------------------------------------------------
Adding new methods to any class is as simple as opening the file containing the class, and adding the method definition.  For instance, say you want the ability to quickly find the average pore size.  You could make a method called Rp_ave() and locate it in GenericGeometry as follows:

.. code-block:: python

   def Rp_ave(self,net):
       return sp.mean(pn.pore_properties['diameter'])

This method will now be available to the rest of the code, or from the command line, as:

.. code-block:: python

   OpenPNM.Geometry.GenericGeometry().Rp_ave(pn)
   
Because theis method was added to the GenericGeometry class it would be available to all geometries by inheritance.  

.. note::
   Of course, this is more typing than simply calculating the average explicity.  It is possible in Python to assign this method to it's own object, which can be accomplished and used as follows:

   .. code-block:: python

      RpAve = OpenPNM.Geometry.GenericGeometry().Rp_ave
      RpAve(pn)

   The second line would return the average pore size.  

===============================================================================
Adding a New Geometry
===============================================================================
Adding a new geometry requires the obvious step of writing the necessary procedures and equations, but it also requires a number of administrative type alterations to the code that allow the new geometry class to register with the rest of the code.

Let's look at the first portion of this task.  A pore network's geometry is defined by the arrangment of pores in space, and by how they are connected by throats.  Although the GenericGeometry class has methods defined for this purpose, these are not implimented; they *must* be implimented in each individual Geometry class.  All of the Geometry classes included with OpenPNM each have their own unique means of defining pores and throats.  





===============================================================================
Manipulating Geometry
===============================================================================

-------------------------------------------------------------------------------
Translate, Scale and Rotate Network
-------------------------------------------------------------------------------
The default geometry generation scheme orients the network relative to [x,y,z] = [0,0,0].  If for any reason the network coordinates must be altered, the GenericGeometry class has several useful tools.

.. automethod:: OpenPNM.Geometry.GenericGeometry.translate_coordinates

.. automethod:: OpenPNM.Geometry.GenericGeometry.scale_coordinates

-------------------------------------------------------------------------------
Stitch Networks
-------------------------------------------------------------------------------
There are several situations where joining or stitching two networks to make a single network is convenient.  One particularly important situation is adding boundary pores to a network.  Given the existence of a cubic network, pn1, of size [10,10,10], boundary pores can be added to a face by generating a second network in memory, pn2, of size [10,10,1].  The new network, pn2, is basically a 2D layer of pores can be added to the face of pn1 to create boundary pores.  Note that both networks have [x,y,z] = [0,0,0] as their origin, so they overlap.  Before peforming the stitch, pn2 should be translated and rotated.  For instance, to attach boundary pores to the x=0 face, the following series of commands would be required:

.. code-block:: python

   OpenPNM.Geometry.GenericGeometry.translate_coords(pn2,[-1,0,0])
   OpenPNM.Geometry.GenericGeometry.stitch(pn1,pn2)
   
This would append the pore properties of pn2 to those of pn1, theyby enlarging pn1.  The pn2 network would remain in memory for subsequent reuse.  

.. automethod:: OpenPNM.Geometry.GenericGeometry.stitch

-------------------------------------------------------------------------------
Extract Sub-Network
-------------------------------------------------------------------------------
...


















