.. _boundary_conditions_example:

===============================================================================
Diffusion Simulations with Various Boundary Conditions
===============================================================================
This example outlines how to perform a variety of transport simulations on pore networks.  It uses the Fickian diffusion algorithm as the basis of the example and demonstrates a number of different boundary conditions specifications.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Generating Network, Geometry, Fluid and Physics
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Start by generating a basic cubic network and the other required components:

.. code-block:: python

    import OpenPNM
    pn = OpenPNM.Network.Cubic(name='test').generate(lattice_spacing=[0.0001],divisions=[10,10,10],add_boundaries=True)
    geo = OpenPNM.Geometry.Stick_and_Ball(network=pn,name='basic')
    geo.regenerate()
    air = OpenPNM.Fluids.Air(network=pn)
    air.regenerate()
    phys = OpenPNM.Physics.GenericPhysics(network=pn,fluid=air,geometry=geo,name='phys')
    phys.add_method(prop='diffusive_conductance',model='bulk_diffusion')
    phys.regenerate()

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Generate an Algorithm Object
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
All algorithms in OpenPNM are independent objects.  The Fickian diffusion algorithm is instantiated as follows:

.. code-block:: python

	alg = OpenPNM.Algorithms.FickianDiffusion(network=pn, name='alg')

	
-------------------------------------------------------------------------------
Apply Dirichlet Conditions to Two Faces
-------------------------------------------------------------------------------

Now this algorithm needs to know something about the boundary conditions which are to be applied.  Let's start by defining Dirichlet conditions on two opposite faces.  This is done by first finding the pores indices that correspond to the two faces.  The generation of cubic networks automatically adds pores to the network with the label 'boundary' as well as a label describing which face they are on.  Let's use 'top' and 'bottom':

.. code-block:: python

	temp_a = pn.get_pore_indices(labels=['top','boundary'],mode='intersection')
	alg.set_pore_info(label='Dirichlet',locations=temp_a)
	temp_b = pn.get_pore_indices(labels=['bottom','boundary'],mode='intersection')
	alg.set_pore_info(label='Dirichlet',locations=temp_b)

The above code adds the label 'Dirichlet' to both 'top' and 'bottom' boundary pores.  The Fickian algorithm looks for this specific label when analyzing and setting up the problem.  Note that the the above code uses the *setter* method associated with the Algorithm object, not the pore network object.  This means that the pore labels will only be applied to this specific algorithm. This is designed to allow multiple algorithms to exist simultaneously without interfering with each other.  The next step is to apply a numerical value to these locations:

.. code-block:: python

	alg.set_pore_data(prop='BCval', data=0.6, locations=temp_a)
	alg.set_pore_data(prop='BCval', data=0.2, locations=temp_b)
	
Note again that the *setter* method of the algorithm was used to keep these boundary conditions isolated from other Algorithms.  The label 'BCval' is a reserved label that that Fickian algorithm, and transport algorithms in general, looks for during the setup.  Once the boundary conditions are specified, the algorithm can be run quite simply as:

.. code-block:: python

	alg.run(active_fluid=air)
	
This runs the algorithm and the results are stored on the Algorithm object.  This is done to prevent simultaneous objects from interfering with each other.  If and when the results of an Algorithm are required by the network model they must be explicitly sent *out* using:

.. code-block:: python

	alg.update()
	
Each Algorithm must subclass the `update()` method so that it sends the correct information out the network and/or fluid.  In the case of the Fickian Algorithm, the 'mole_fraction' of the active_fluid is stored on the Fluid object in question.  Running a different version of the Algorithm and calling `update()` will overwrite any previous values.  The results of this simulation should produce the following visualization (done in Paraview):

.. image:: BC1.png
	
-------------------------------------------------------------------------------
Apply Neumann Conditions to a Group of Internal Pores
-------------------------------------------------------------------------------

The code below sets the total rate leaving a group of pores cumulatively.  Note that the same Algorithm object is used (`alg`), so the Dirichlet boundary conditions applied in the previous step still exist.  The lines below define a group of 5 pores which are generating mass at a set rate, which is accomplished by creating a 'Neumann_rate_group' label and placing the numerical value of the rate in 'BCval' array.  

.. code-block:: python

	temp_c = [500,501,502,503,504]
	alg.set_pore_info(label='Neumann_rate_group',locations=temp_c)
	alg.set_pore_data(prop='BCval',data=5e-7,locations=temp_c)
	alg.run(active_fluid='air')
	alg.update()

This results in the image below, where a region of high concentration can be seen in the core of the domain due to the mass production: 

.. image:: BC2.png

-------------------------------------------------------------------------------
Apply Neumann Conditions in Several Pores Individually
-------------------------------------------------------------------------------

One of the options for specifying Neumann conditions is to apply the same rate to multiple pores.  Begin by removing some of the conditions applied above, then set a few pores on the 'bottom' face to each have the same specific rate.

.. code-block:: python

	alg.set_pore_info(label='Neumann_rate_group',locations=temp_c,mode='remove')  # This removes label from pores
	alg.set_pore_info(label='Dirichlet',locations=temp_b,mode='remove')
	alg.set_pore_info(label='Neumann_rate_single',locations=temp_b)
	alg.set_pore_data(label='BCval',data=1e-10,locations=temp_b)
	alg.run(active_fluid='air')
	alg.update()
	
This results in image below.  Notice that the concentration on the inlet face is not uniform, and that the smaller pores have a somewhat higher concentration (darker red), which is necessary if their flux is the be the same as larger, more conductive pores.

.. image:: BC3.png







