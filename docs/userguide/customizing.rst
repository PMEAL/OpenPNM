.. _customizing:

================================================================================
Customizing
================================================================================

.. contents:: Page Contents
    :depth: 3

The OpenPNM framework was designed specifically with extensibility and customization in mind.  Every user will apply OpenPNM to a unique problem, and will therefore require unique pore scale models, phase properties, algorithms and so on.

There are two ways to customize OpenPNM.  The first is to download the source code and *hack* it.  With this approach it is possible to create your own sub-classes, add pore-scale models, define new topology generators, and to add or change OpenPNM methods.  The other approach is to install OpenPNM in your Python PATH (using ``pip install openpnm``) as with any other package such as Scipy, and write custom functions in a separate 'working directory' rather than in the source code.  In this scenario, you can perform all of the same customizations as the first approach, with the exception of changing OpenPNM's native methods (in fact even this is possible, but that's another story).  The second approach is the recommended way for several reasons.  It avoids accidental changes to the framework, it allows users to keep their 'projects' compartmentalized, and it is much easier for users to contribute their work to OpenPNM project (which we highly encourage) since sections can be 'cut and pasted' or 'merged' into the framework cleanly.  The second approach will be explained in detail below.

The following discussions assume that all custom files will be stored in a folder called ``my_pnm``, that will be the 'working directory'.

-------------------------------------------------------------------------------
Creating Custom Models
-------------------------------------------------------------------------------
In the working directory, place a file called 'my_models.py' (or any name you wish).  This file will be the home of all the custom models that will be created. Models that are in this file can be added to *any* object using the ``add_model`` command.  The *'models'* mechanism in OpenPNM was designed to be as straight-forward as possible, for users without any experience in object oriented coding.  Each model is simply a 'function' definition, with no inheritance, classes or any unnecessary complications.

Let's create a model called 'surface_roughness' that calculates the surface area of a pore accounting for surface roughness, that accepts a single 'roughness parameter' and is a function of pore size.  Start by writing a function in the 'my_models.py' file that simply accepts the 'roughness_parameter' and simply returns it:

.. code-block:: python

    def surface_roughness(roughness_parameter,**kwargs):
    	return roughness_parameter

The use of the 'kwargs' argument is essential.  Without diving into the details, 'kwargs' will collect all arguments that are sent to function but not used.  Without this, the function will break since the ``add_model`` method sends many arguments 'just-in-case'.

Now, we can see this model in action by creating a script in the working directory as follows:

.. code-block:: python

	import my_models
	a = my_models.surface_roughness(roughness_parameter=2.2)
	print(a)
	2.2

The next step is to have this model to calculate something useful.  Assume that surface area due to roughness scales with the projected or smooth surface area as some function, which means this function will need access to the 'pore.diameter' information, which is stored on Geometry objects.  Thus the relevant Geometry object must be sent as an argument:

.. code-block:: python

	def surface_roughness(geometry,roughness_parameter,**kwargs):
		P_diam = geometry['pore.diameter']
		projected_area = 4*3.14159*(P_diam/2)**2
		rough_area = projected_area**roughness_parameter
		return rough_area

It was noted above that the ``add_model`` method sent in several extra arguments 'just-in-case'.  Among these arguments are the object from which the method is called.  Since 'surface_roughness' will be attached to a Geometry object, this function will receive it as 'geometry' automatically.  We can now update our script:

.. code-block:: python

	import OpenPNM
	import scipy as sp
	import my_models

	#Generate a simple Cubic Network
	pn = OpenPNM.Network.Cubic(shape=[3,3,3])

	#Generate an 'empty' Geometry with no properties
	geom = OpenPNM.Geometry.GenericGeometry(network=pn,pores=pn.pores(),throats=pn.throats())

	#Assign random pores diameters between 0 and 40
	geom['pore.diameter'] = sp.rand(pn.Np,)*40

	#Assign model to geometry
	geom.add_model(propname='pore.surface_area',
					model=my_models.surface_roughness,
					roughness_parameter=2.2)

The print-out of 'geom 'reveals that indeed the model has been added:

print(geom)
------------------------------------------------------------
OpenPNM.Geometry.GenericGeometry: 	GenericGeometry_4rhgW
------------------------------------------------------------
#     Properties                          Valid Values
------------------------------------------------------------
1     pore.diameter                          27 / 27
2     pore.surface_area                      27 / 27
------------------------------------------------------------
#     Labels                              Assigned Locations
------------------------------------------------------------
1     pore.all                            27
2     throat.all                          54
------------------------------------------------------------

The same approach can be used to create models for pore-scale Physics or for calculating fluid properties that are not included with OpenPNM.

-------------------------------------------------------------------------------
Using non-Default Property Names
-------------------------------------------------------------------------------
In the 'surface_roughness' example above, the function assumed that pore diameter data would be found under the 'pore.diameter' dictionary key.  If for some reason, there were multiple different definitions of 'pore diameter', then they might be stored as 'pore.diameter_inscribed', and 'pore.diameter_hydraulic', etc.  To allow the 'surface_roughness' function to be applied to any arbitrary pore diameter, it should be rewritten as:

.. code-block:: python

	def surface_roughness(geometry, roughness_parameter, pore_diameter='pore.diameter', **kwargs):
		P_diam = geometry[pore_diameter]
		projected_area = 4*3.14159*(P_diam/2)**2
		rough_area = projected_area**roughness_parameter
		return rough_area

Note that *pore_diameter* is now an argument name, which defaults to 'pore.diameter'.  Different *pore diameters* can be specified when calling ``add_model``:

.. code-block:: python

	#Assign model to geometry
	geom.add_model(propname='pore.surface_area',
				   model=my_models.surface_roughness,
				   pore_diameter = 'pore.diameter_inscribed',
				   roughness_parameter=2.2)

All of the models provide with OpenPNM allow for this sort of non-default argument names, and it will make your custom models more general if you follow this practice.

-------------------------------------------------------------------------------
Creating a Customized Subclass
-------------------------------------------------------------------------------
Another important way to customize OpenPNM is to create custom subclasses of the various objects.  For instance, OpenPNM comes with a few basic Geometry subclasses that return pore-scale geometric properties representative of various materials, or common fluids.  Creating a custom subclass is only slightly more complicated than writing custom models.

Let's create a Geometry subclass that is representative of Berea Sandstone.  Start by creating a file in the 'working directly' (assume it's called 'my_geometries').  In this file we need define our 'class', which will inherit from OpenPNM.Geometry.GenericGeometry:

.. code-block:: python

	import OpenPNM
	class BereaSandstone(OpenPNM.Geometry.GenericGeometry):
		def __init__(self,**kwargs):
			super(berea_sandstone,self).__init__(**kwargs)

The above is a basic template that is no different than GenericGeometry yet.  The important thing to notice here is the the ``__init__`` of the parent class is invoked using the ``super`` method.  This means that all arguments passed to ``BereaSandstone`` are bundled into *'kwargs'* and passed to **GenericGeometry**, which will run all of the tasks that are necessary for OpenPNM objects to work, such as registering this custom Geometry with the Network.

The next step is to actually customize the class.  In OpenPNM, all the subclasses of Geometry, Phase and Physics are literally just a collection of 'models' with appropriate parameters to reproduce a specific material, fluid or set of physics.  The BereaSandstone class then just needs a set of suitable 'models':

.. code-block:: python

	import OpenPNM
	class BereaSandstone(OpenPNM.Geometry.GenericGeometry):
		def __init__(self,**kwargs):
			super(berea_sandstone,self).__init__(**kwargs)

		mod = OpenPNM.Geometry.models.pore_misc.random
		self.add_model(propname='pore.seed',
						model=mod)

		mod = OpenPNM.Geometry.models.pore_diameter.sphere
		self.add_model(propname='pore.diameter',
						model=mod,
						psd_name='weibull_min',
						psd_shape=2.5,
						psd_loc=4e-4,
						psd_scale=4e-4)


The first of the above two models creates a property called 'pore.seed', which is just a list of random numbers that will be used to seed the pore size distribution.  The second model uses the Scipy.stats package to generate 'pore.diameter' values from the 'weibull_min' distribution using the given parameters.

--------------------------------------------------------------------------------
Creating Customized Networks
--------------------------------------------------------------------------------

Unlike Geometry, Phase and Physics objects, a Network object requires more than a collection of calls to ``add_model``.  The Network object must provide the 'pore.coords' and 'throat.conns' properties.  The 'pore.coords' is fairly straightforward, as it's just an Np x 3 list of [x,y,z] coordinates for each pore in the Network.  The 'throat.conns' list is much more difficult to produce.  This list is an Nt x 2 list of pairs of connected pore, such as [P1,P2].  OpenPNM comes with two main Network classes: Cubic and Delaunay.  The Cubic class connects each pore to it's immediate 6 neighbors on a cubic lattice, while the Delaunay class places pores randomly in space and determines connections via a Delaunay tessellation.  There are endless possible topology generation schemes that one may wish to develop.

The approach used to subclass GenericGeometry above would also work for Networks, but there is one additional consideration.  Every object must have a 'pore.all' and a 'throat.all' array so that they function properly.  The Network generation must therefore, produce these two arrays as well as the 'pore.coords' and 'throat.conns' described above.

--------------------------------------------------------------------------------
Creating Customized Algorithms
--------------------------------------------------------------------------------

Algorithms can also be customized as described above.  The GenericAlgorithm has a few additional methods that are meant to be implemented by subclasses, such as ``return_results``.  The intention of this method is to send the pertinent results of a calculation 'out' of the Algorithm object and to the correct object in the simulation.  This step is handy, but is not actually necessary.  One can of course manually transfer data from an Algorithm to a Phase, for instance with:
