
.. _cubic-example:

===============================================================================
Generate and Use a Cubic Network
===============================================================================

-------------------------------------------------------------------------------
Generating Geometry
-------------------------------------------------------------------------------
Start by importing OpenPNM

.. code-block:: python

    import OpenPNM
    
Next, setup a dictionary containing all the desired network parameters. 

.. code-block:: python

    params = {
    'domain_size'           : [0.001,0.001,0.0004],  #physical network size [meters]
    'divisions'             : [], #Number of pores in each direction
    'lattice_spacing'       : [.0001],  #spacing between pores [meters]
    'num_pores'             : 1000, #This is used for random networks where spacing is irrelevant
    'template'              : template, #This is used for the Template based network generation
    'psd_info'   : {'name'  : 'weibull_min', #Each statistical package takes different params, so send as dict
                    'shape' : 1.5,
                    'loc'   : 6e-6,
                    'scale' : 2e-5},
    'tsd_info'   : {'name'  : 'weibull_min',
                    'shape' : 1.5,
                    'loc'   : 6e-6,
                    'scale' : 2e-5},
    'btype'                 : [0,0,0],  #boundary type to apply to opposing faces [x,y,z] (1=periodic)
    }

Once the parameters are specified, a network called 'pn' can be generated with:

.. code-block:: python

    pn = OpenPNM.Geometry.Cubic(loglevel=10).generate(**params)

Note that Cubic() accepts arguments that control the overall behavior of the class, such as setting the loglevel (for very verbose logging output set the loglevel=10, and for only urgent warnings and errors set loglevel=50).

Now a cubic network of 1000 pores is stored in the object named 'pn'.  To see what sort of properties have been generated for this network, you can ask it to print an overview of itself as follows:

.. code-block:: python

    pn.print_overview()
    
This indicates all the pore and throat properties that are currently associated with the network.  You'll notice that they are exclusively geometric in nature, which is normal since no simulations have yet been run.  

-------------------------------------------------------------------------------
Using Pore Scale Physics
-------------------------------------------------------------------------------
To perform a capillary pressure curve simulation we must first generate a throat property that describes the capillary pressure required for the non-wetting fluid to invade as a function of size, wettability, fluid properties and so on.  OpenPNM comes with a Physics Module which contains numerous models and equations for calculating such properties.  We can generate the most basic estimate of capillary entry pressure by using the Washburn equation for cylindrical tubes:

.. code-blocK:: python

	OpenPNM.Physics.CapillaryPressure.Washburn(pn, sigma = 0.72, theta = 120)
	
This method calculates the capillary entry pressure of each throat in the network based on the surface tension and wettability information provided.  The results of this calculation are stored in throat_conditions['Pc_entry'].  Note that this this is a 'condition' rather than a 'property' because it depends on which fluid is used so it's not an intrinsic network property.  Also, note that this 'condition' was stored in :code:`throat_conditions` since capillary invasion is controlled by throats not pores.  


-------------------------------------------------------------------------------
Running Simulations
-------------------------------------------------------------------------------






-------------------------------------------------------------------------------
Visualizing with Paraview
-------------------------------------------------------------------------------







