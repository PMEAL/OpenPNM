###############################################################################
Getting Started
###############################################################################

===============================================================================
Example: Working with Scripts
===============================================================================
1.  Open Spyder
2.  In Spyder, select the 'Browse a Working Directory' button and locate your OpenPNM folder and select.  Also click the 'Set as console's working directory' button.
3.  Open a new blank .py file in the editor window if there is not one already.  

The first thing you must do is tell Spyder to import the OpenPNM code so you have access to the functions and methods, so add

.. code-block:: python

   >>> import OpenPNM

This imports the OpenPNM package with all it's algorithms and data storage functions.

Next, you'll want to generate a network.  This is accomplished by first creating a generator, which we'll call gn1:

.. code-block:: python
   
   >>> gn1 = OpenPNM.Geometry.Cubic()
   
Now a generator has been initialized (with default parameters) and is ready to start spawning networks.  This can be achieved by:

.. code-block:: python
   
   >>> pn1 = gn1.generate()

To view the properties of pn1, use the print command:

.. code-block:: python

    >>> print pn1
    ==================================================
    Overview of network properties
    --------------------------------------------------
    Basic properties of the network
    - Number of pores:   81
    - Number of throats: 108
    
    Pore properties:
    	diameter            float64             (81,)               
    	numbering           int32               (81,)               
    	volume              float64             (81,)               
    	seed                float64             (81,)               
    	coords              float64             (81, 3)             
    	type                int8                (81,)               
    Throat properties:
    	volume              float64             (108,)              
    	diameter            float64             (108,)              
    	numbering           int32               (108,)              
    	connections         int32               (108, 2)            
    	length              float64             (108,)              
    	seed                float64             (108,)              
    	type                int8                (108,)  

To create a second network, you must initiate a second geometry object:

.. code-block:: python
   
   >>> pn2 = OpenPNM.Geometry.Cubic().generate()
   
This network will be identical in most asepcts but will have different pore and throat sizes due the different random seeds used.  This functionality is useful for running simulations on multiple realizations for the same network for statistical significance.  

===============================================================================
Example: Querying the Network
===============================================================================




