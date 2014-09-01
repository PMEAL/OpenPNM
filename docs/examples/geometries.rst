.. _geometries_example:

===============================================================================
Generating Geometries
===============================================================================

This example demonstrates how to create and manipulate *Geometry* objects, including applying multiple Geometries to the network, removing geometries, reassigning geometries to different locations, and so on. 

Start by importing OpenPNM and creating a Network object:

.. code-block:: python

	import OpenPNM
	pn = OpenPNM.Network.Cubic.empty(name='net',loglevel=20,dims=[10,10,10])
	pn.add_boundaries()
	
.. note:: The Importance of Labels
	
	Labels play a central role in OpenPNM.  It allows for pores and throats to be *tagged*, or *categorized* into certain groups, which makes it very convenient to select a set of pores or throats as needed.  For instance, the Cubic network generation adds several labels to the Network by default, such as 'top', 'front', etc.  To select all the pores on the 'top' of the network, it is simply a matter of finding where the label of 'top' has been applied.  Of course, OpenPNM incorporates several tools for this.  To see a list (printed to the command lines) of all the labels currently applied to the network use ``pn.labels()``.  The process of selecting pores with a certain label (or labels) is demonstrated below:

	.. code-block:: python
	
		pores = pn.pores(labels='top')  # Finds top face
		pores = pn.throats(labels=['top','front'], mode='intersection')  # Finds top-front edge
	
	A similar approach is used to find throats.  More details on the usage these functions and their options can be found in their documentation.  


 
 Also, since the ``pn.add_boundaries()`` function was called after instantiating the ``pn`` object,  the domain is padded with a layer of boundary pores that are labelled as 'boundary'. 