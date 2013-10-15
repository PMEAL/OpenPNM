
.. _cubic-example:

===============================================================================
Generate Cubic Network and Simulate a Drainage Capillary Pressure Curve
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
    'domain_size'           : [10,10,10],  #physical network size [meters]
    'divisions'             : [], #Number of pores in each direction
    'lattice_spacing'       : [1],  #spacing between pores [meters]
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
In order to run a drainage simulation it is necessary to specify the inlet sites from which the invasion of non-wetting fluid proceeds.  There are several possibilities here depending on what sort experiment is being simulated.  For mercury intrusion porosimetry (MIP), the non-wetting fluid invades the sample from all sides.  During the `generate()` step, boundary pores were added on all sides of the network and give a `'type'` value > 0 (0 indicates an internal pores).  To specify invasion from all faces, the inlets can be set to all boundary pores:

.. code-block:: python

	mask = pn.pore_properties['type']>0
	inlets = pn.pore_properties['numbering'][mask]

The simulation can be run as:

.. code-block:: python

	OpenPNM.Algorithms.OrdinaryPercolation(loglevel = 10).run(net = pn, npts = 50, inv_sites = inlets)
	
The resulting capillary pressure curve can be visualized by sending the network (pn) to the custom built plot command available in the Visualization module:

.. code-block:: python

	OpenPNM.Visualization.Plots.Capillary_Pressure_Curve(pn)

The capillary pressure curve should like something like this:

.. plot::
	
	import matplotlib.pyplot as plt
	import OpenPNM
	pn = OpenPNM.Geometry.Cubic(loglevel=10).generate()
	OpenPNM.Physics.CapillaryPressure.Washburn(pn, sigma = 0.72, theta = 120)
	mask = pn.pore_properties['type']>0
	inlets = pn.pore_properties['numbering'][mask]
	OpenPNM.Algorithms.OrdinaryPercolation(loglevel = 10).run(net = pn, npts = 50, inv_sites = inlets)
	plt.hist(pn.pore_conditions['Pc_invaded'])
   

-------------------------------------------------------------------------------
Visualizing with Paraview
-------------------------------------------------------------------------------







