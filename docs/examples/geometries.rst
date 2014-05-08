.. _geometries_example:

===============================================================================
Generating Geometries
===============================================================================

This example demonstrates how to create and manipulate *Geometry* objects, including applying multiple Geometries to the network, removing geometries, reassigning geometries to different locations, and so on. 

Start by importing OpenPNM and creating a Network object:

.. code-block:: python

    import OpenPNM
    pn = OpenPNM.Network.Cubic(name='test').generate(lattice_spacing=[0.0001],divisions=[10,10,10],add_boundaries=True)
	
.. note:: The Importance of Labels
	
	Labels play a central role in OpenPNM.  It allows for pores and throats to be *tagged*, or *categorized* into certain groups, which makes it very convenient to select a set of pores or throats as needed.  For instance, the Cubic network generation add several labels to the Network by default, such as 'top', 'front', etc.  To select all the pores on the 'top' of the network, it is simply a matter of finding where the label of 'top' has been applied.  Of course, OpenPNM incorporates several tools for this.  To see a list (printed to the command lines) of all the labels currently applied to the network use ``list_pore_labels`` and ``list_throat_labels``.  The process of selecting pores with a certain label (or labels) is demonstrated below:

	.. code-block:: python
	
		pores = pn.get_pore_indices(labels='top')  # Finds top face
		pores = pn.get_pore_indices(labels=['top','front'], mode='intersection')  # Finds top-front edge
	
	A similar approach is used to find throats.  More details on the usage these functions and their options can be found in their documentation.  


 
 Also, since the add_boundaries flag was set to True, the ``generate`` function will have padded the domain with a layer of boundary pores and labeled them as 'boundary'. 