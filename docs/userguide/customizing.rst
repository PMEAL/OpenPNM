.. _customizing:

================================================================================
Customizing
================================================================================

.. contents:: Page Contents
    :depth: 3

The OpenPNM framework was designed with extensibility in mind.  Every user will apply OpenPNM to a unique problem, and will therefore require unique pore scale models, phase properties, algorithms and so on.

There are two ways to customize OpenPNM.  The first is to download the source code and *hack* it.  With this approach it is possible to create your own sub-classes, add pore-scale models, define new topology generators, and to add or change OpenPNM methods.  The other approach is to install OpenPNM in your Python PATH (using ``pip install openpnm``) as with any other package such as Scipy, and write custom functions in a separate 'working directory' rather than in the source code.  The second approach is the recommended way for several reasons.  It avoids accidental changes to the framework, it allows users to keep their 'projects' compartmentalized, and it is much easier for users to contribute their work to OpenPNM project since sections can be merged into the framework cleanly.  The second approach will be explained in detail below.

The following discussions assume that all custom files will be stored in a folder called ``'my_pnm'``, that will be the 'working directory'.

-------------------------------------------------------------------------------
Creating Custom Models
-------------------------------------------------------------------------------
In the working directory, place a file called 'my_models.py'.  This file will be the home of all the custom models that will be created. Models that are in this file can be added to *any* object using the ``add_model`` command.  The *'models'* mechanism in OpenPNM was designed to be as straight-forward as possible, so each model is simply a function definition.

Let's create a model called 'surface_roughness' that calculates the surface area of a pore accounting for surface roughness, that accepts a single 'roughness parameter' and is a function of pore size.  Start by writing a function in the 'my_models.py' file that simply accepts the 'roughness_parameter' and returns it:

.. code-block:: python

  def surface_roughness(roughness_parameter):
    return roughness_parameter

Now, we can see this model in action by creating a script in the working directory as follows:

.. code-block:: python

	import my_models
	a = my_models.surface_roughness(roughness_parameter=2.2)
	print(a)
	2.2

The next step is to have this model calculate something useful.  Assume that surface area due to roughness scales with the projected or smooth surface area as some function, which means this function will need access to the 'pore.diameter' information, which is stored on Geometry objects.  Thus the relevant Geometry object must be sent as an argument.  OpenPNM assumes object for which the data is being calculated is passed in as `target`, which makes all function calls general:

.. code-block:: python

  def surface_roughness(target, roughness_parameter):
  	P_diam = target['pore.diameter']
  	projected_area = 4*3.14159*(P_diam/2)**2
  	rough_area = projected_area**roughness_parameter
  	return rough_area

We can now update our script:

.. code-block:: python

	import openpnm as op
	import scipy as sp
	import my_models

	#Generate a simple Cubic Network
	pn = op.network.Cubic(shape=[3,3,3])

	#Generate an 'empty' Geometry with no properties
	geom = op.geometry.GenericGeometry(network=pn,pores=pn.pores(),throats=pn.throats())

	#Assign random pores diameters between 0 and 40
	geom['pore.diameter'] = sp.rand(pn.Np,)*40

	#Assign model to geometry
	geom.add_model(propname='pore.surface_area',
					       model=my_models.surface_roughness,
					       roughness_parameter=2.2)

The same approach can be used to create models for pore-scale Physics or for calculating fluid properties that are not included with OpenPNM.

-------------------------------------------------------------------------------
Using non-Default Property Names
-------------------------------------------------------------------------------
In the 'surface_roughness' example above, the function assumed that pore diameter data would be found under the 'pore.diameter' dictionary key.  If for some reason, there were multiple different definitions of 'pore diameter', then they might be stored as 'pore.diameter_inscribed', and 'pore.diameter_hydraulic', etc.  To allow the 'surface_roughness' function to be applied to any arbitrary pore diameter, it should be rewritten as:

.. code-block:: python

  def surface_roughness(target, roughness_parameter, pore_diameter='pore.diameter'):
	   P_diam = target[pore_diameter]
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
Creating a Basic Subclass: Predefining Pore-Scale Models
-------------------------------------------------------------------------------
Another way to customize OpenPNM is to create custom subclasses.  This approach is best if you wish to apply several pore-scale models to a single object, since it let's you collect them all into one place.  This is done by starting with one of OpenPNMs existing Generic classes and adding a suite of pore-scale model definitions in the `init` stage.  This is how OpenPNM classes such as Water, and StickAndBall operate.  They do not actually overload or add any method, but just act as a collection of pre-set pore-scale models.

For example, let's create a custom Phase object for an oil with temperature dependent viscosity and density.  The following class definition can be added to ``'my_classes.py'`` in the same directory as ``'my_models.py'``.

.. code-block:: python

  from openpnm.phases import GenericPhase


  class Oil(GenericPhase):
      def __init__(self, **kwargs):
          super().__init__(**kwargs)  # This is a python thing, and calls the init of the parent class

    		  self.add_model(propname='pore.viscosity',
             						 model=op.models.misc.polynomial,
                         a=[10000, -111, 2],
                         prop='pore.temperature')

    			self.add_model(propname='pore.density',
    			 			         model=op.models.misc.linear,
    						         prop='pore.temperature',
                         m=2200, b=-20)

The first of the above two models creates a property called 'pore.viscosity', which using a polynomial function to describe the dependence on temperature (indicated by the ``prop`` argument).  The second model is similar but using a linear fitting to describe the density.  The values used for the coefficients of these two models will dicatate the final physical properties of the Phase, so this is now a custom phase.

Thus when you change the values of 'pore.temperature' on an instance of ``Oil``, then call `regenerate_models`, these two models will be run and will look at the current value in 'pore.temperture'.

--------------------------------------------------------------------------------
Creating Customized Networks
--------------------------------------------------------------------------------

Unlike Geometry, Phase and Physics objects, a Network object requires more than a collection models.  In fact, Networks typically have no models.  Creating a custom Network type is all about defining the locations of the pores, and the connections between the throats.  Consider the following basic graph:

::

    4 ----- 3
    | \   / |
    5 - 0   |
    |     \ |
    2 - 6 - 1

Implementing this as a 'custom' Network can be done as by noting that the pore coordinates are:

.. code-block:: python

    coords = [[1, 1, 0],
              [2, 0, 0],
              [0, 0, 0],
              [2, 2, 0],
              [0, 2, 0],
              [0, 1, 0],
              [1, 0, 0]]

And the connections are:

.. code-block:: python

    conns = [[0, 1],
             [0, 3],
             [0, 4],
             [0, 5],
             [1, 3],
             [1, 6],
             [2, 5],
             [2, 6],
             [3, 4],
             [4, 5]]

A custom Network can be created by passing these two arrays in to a GenericNetwork class:

.. code-block:: python

    pn = op.network.GenericNetwork(coords=coords, conns=conns)

The GenericNetwork class is able to handle topologies of any sort, so you can be as imaginitive as you wish when defining the network.  There are several rules and assumptions about how the ``conns`` and ``coords`` data must be formatted, which are descibed in :ref:`topology`

To create an actual Network subclass, you can add the following to your ``'my_networks.py'`` file:

.. code-block:: python

  from openpnm.network import GenericNetwork


  class MyNetwork(GenericNetwork):
      def __init__(self, **kwargs):
          # Place pore locations and throat connecting generation code here
          super().__init__(conns=conns, coords=coords, **kwargs)


--------------------------------------------------------------------------------
Creating Customized Algorithms
--------------------------------------------------------------------------------

Algorithms can also be customized as described above.  The GenericAlgorithm has a few additional methods that are meant to be implemented by subclasses, such as `return` and `reset`.  The intention of this method is to send the pertinent results of a calculation 'out' of the Algorithm object and to the correct object in the simulation.  This step is handy, but is not actually necessary.  One can of course manually transfer data from an Algorithm to a Phase, for instance with:
