
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
Creating Fluids
-------------------------------------------------------------------------------
.. code-block::python
	air_recipe = {       'name': 'air',
						   'Pc': 3.771e6, #Pa
						   'Tc': 132.65,  #K
						   'MW': 0.0291,  #kg/mol
				  'diffusivity': {'method': 'Fuller',
									  'MA': 0.03199,
									  'MB': 0.0291,
									  'vA': 16.3,
									  'vB': 19.7},
					'viscosity': {'method': 'Reynolds',
									  'uo': 0.001,
									   'b': 0.1},
				'molar_density': {'method': 'ideal_gas',
									   'R': 8.314},
			  'surface_tension': {'method': 'constant',
								   'value': 0},
				'contact_angle': {'method': 'na'},
	}
	water_recipe = {     'name': 'water',
						   'Pc': 2.206e6, #Pa
						   'Tc': 647,     #K
						   'MW': 0.0181,  #kg/mol
				  'diffusivity': {'method': 'constant',
								   'value': 1e-12},
					'viscosity': {'method': 'constant',
								   'value': 0.001},
				'molar_density': {'method': 'constant',
								   'value': 44445},
			  'surface_tension': {'method': 'Eotvos',
									   'k': 2.25e-4},
				'contact_angle': {'method': 'constant',
								   'value': 120},
	}

Now that the fluids *recipes* are defined they can be passed to the `create()` method of the fluids module:

.. code-block::python
	#Create fluids
	air = OpenPNM.Fluids.GenericFluid(loglevel=50).create(air_recipe)
	water= OpenPNM.Fluids.GenericFluid(loglevel=50).create(water_recipe)
	#Set desired base conditions in the Fluids
	air.pore_conditions['temperature'] = 353
	air.pore_conditions['pressure'] = 101325
	water.pore_conditions['temperature'] = 353
	water.pore_conditions['pressure'] = 101325
	#Update Fluids to the new conditions
	water.regenerate()
	air.regenerate()

-------------------------------------------------------------------------------
Using Pore Scale Physics
-------------------------------------------------------------------------------
To perform a capillary pressure curve simulation we must first generate a throat property that describes the capillary pressure required for the non-wetting fluid to invade as a function of size, wettability, fluid properties and so on.  OpenPNM comes with a Physics Module which contains numerous models and equations for calculating such properties.  We can generate the most basic estimate of capillary entry pressure by using the Washburn equation for cylindrical tubes:

.. code-blocK:: python

	OpenPNM.Physics.CapillaryPressure.Washburn(pn, water)
	
This method calculates the capillary entry pressure of each throat in the network based on the surface tension and wettability information provided.  The results of this calculation are stored as part of the water fluid undet throat_conditions['Pc_entry'].  Note that this 'condition' was stored in :code:`throat_conditions` since capillary invasion is controlled by throats not pores.  

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

	OpenPNM.Visualization.Plots.Capillary_Pressure_Curve(water)

   

-------------------------------------------------------------------------------
Visualizing with Paraview
-------------------------------------------------------------------------------







