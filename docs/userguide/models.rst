.. _models:

===============================================================================
Pore Scale Models
===============================================================================
Models are one of the most important aspects of OpenPNM, as they allow the user to specify a 'model' for calculating 'pore.volume', rather than just entering numerical values into a geometry_object['pore.volume'] array.  It is mainly through customized models that users can tailor OpenPNM to a specific situation, though OpenPNM includes a variety of pre-written models.  These are stored under each Module in a folder called 'models'.  For instance, Geometry.models.pore_diameter contains several methods for calculating pore diameters.  For an example on creating custom models see :ref:`Customizing OpenPNM<customizing>`.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Models Dictionary
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Each Core object has a ``models`` attribute where all information about pore-scale models are stored.  Upon instantiation of each ``Core`` object, a ``ModelsDict`` object is stored in its ``models`` attribute.  The ``ModelsDict`` class is a subclass of the Python ``dict`` class, which has several features added for dealing specifically with models.  Each ``model`` and its associated arguments are wrapped in a ``GenericModel`` dictionary, then added to the ``ModelsDict`` under the specified ``propname``.  When a model is run the values is produces are automatically stored in the ``Core`` objects dictionary under the same specified ``propname``.  

Adding a model to an object is done as follows:

(1) A handle to the desired model is retrieved, either from the included OpenPNM model libraries, or from a file containing the users custom models.
(2) The model is attached to the target object using ``add_model``.

This process is demonstrated by added a random pore seed model to a Geometry object:

.. code-block:: python

	geom = OpenPNM.Geometry.GenericGeometry()  # Creates an empty Geometry object
	mod = OpenPNM.Geometry.models.pore_misc.random  # Get a handle to the desired model
	geom.add_model(propname='pore.seed',  # Specify the name of the property output by the model
	               model=mod,             # Assign model to the object
				   seed=0)                # Specify any arguments required by the model
	
The *propname* and *model* arguments are required by the ``add_model`` method, but the *seed* argument is passed on the model, and it specifies the initialization value for the random number generator.  

By default, the ``add_model`` method runs the model and places the data in the Core object's dictionary under the given ``propname``. It also wraps the handle to the model and all the given parameters into a ``GenericModel`` dict (described below), then saves its in the ``ModelsDict`` under the same ``propname``.  

In order to recalculate the data the model stored in the private dictionary must be rerun.  This is accomplished with the ``regenerate`` method.  This method takes an optional list of *propnames* that should be regenerated.  It should also be pointed out that models are regenerated in the order that they were added to the object so some care must be taken to ensure that changes in property values cascade through the object correctly.  The ``ModelsDict`` class has functions for updating the model order.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Generic Model
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++