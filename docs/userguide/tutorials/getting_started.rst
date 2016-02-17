.. _getting_started:

###############################################################################
Tutorial 1 of 3: Getting Started with OpenPNM
###############################################################################

As usual, start by importing the OpenPNM package, and the Scipy package:

>>> import OpenPNM
>>> import scipy as sp

===============================================================================
Building a Cubic Network
===============================================================================

Start by generating a **Network**.  This is accomplished by choosing the desired network topology (e.g. cubic), then calling its respective method in OpenPNM with the desired parameters:

>>> pn = OpenPNM.Network.Cubic(shape=[10, 10, 10], spacing=0.0001)

This generates a topological network and stores it in variable ``pn``.  This network contains pores at the correct spatial positions and connections between the pores according the specified topology (but without boundary pores).  The ``shape`` argument specifies the number of pores in the [X, Y, Z] directions of the cube.  Networks in OpenPNM are alway 3D dimensional, meaning that a 2D or 'flat' network is still 1 layer of pores 'thick' so [X, Y, Z] = [20, 10, 1].  The ``spacing`` argument controls the center-to-center distance between pores.  Although OpenPNM does not currently have a dimensional units system, we *strongly* recommend using SI throughout.

The Network object has numerous methods that can be used to query the topological properties:

>>> pn.num_pores()
1000
>>> pn.num_throats()
2700
>>> pn.find_neighbor_pores(pores=[1])  # Find neighbors of pore 1
array([  0,   2,  11, 101])
>>> pn.find_neighbor_throats(pores=[1, 2])  # Find throats connected to pores 1 and 2
array([   0,    1,    2,  901,  902, 1801, 1802])

There are several more such topological query method available on the object such as ``find_nearby_pores`` and ``find_connecting_throat``.  A full list is given in the detailed documentation <HERE>, along with an explanation of each argument and some helpful examples.

Another important feature is the use of *labels* on pores and throats.  Applying a label to a set of special pores allows for easy retrieval of these pores for later use.  For instance, during the generation of a Cubic network, the faces are automatically labeled.  The following illustrates how to use labels:

>>> pn.labels(pores=[1])  # Find all labels applied to pore 1
['pore.all', 'pore.front', 'pore.internal', 'pore.left']
>>> pn.pores(labels=['front', 'left'], mode='intersection')
array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])

This last command has an argument ``mode``.  In this case it means the set logic to apply to the query (i.e. only return pores that are labeled both 'front' and 'left').

.. note::

	OpenPNM includes extensive help within the code that explains the use of each method, the meaning of each argument and option, and even gives examples.  Within the Spyder IDE these can be viewed in the *object inspector* by typing *ctrl-i* while the cursor is on the method name (in either the editor or the console).  Importantly, Numpy, Scipy, Matplotlib, Pandas and most other scientific packages contain similar detailed documentation.  These are indispensible and even seasoned OpenPNM coders rely on this documentation.

===============================================================================
Initialize and Build a Geometry Object
===============================================================================

The **Network** does not contain any information about pore and throat sizes at this point.  The next step, then, is to create a **Geometry** object to calculate the desired geometrical properties.

>>> geom = OpenPNM.Geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)

This statement contains three arguments: ``network`` tells the **Geometry** object which **Network** it is associated with.  ``pores`` and ``throats`` indicate which locations in the **Network** where this **Geometry** object will apply.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Add Desired Properties to Geometry
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

This freshly instantiated **Geometry** object contains no geometric properties as yet because we chose to use the *GenericGeometry* class.
-------------------------------------------------------------------------------
Direct Assignment of Static Values
-------------------------------------------------------------------------------

Let's start by assiging diameters to each pore from a random distribution, spanning 10 um to 100 um.  The upper limit arises because the ``spacing`` of the *Network* was set to 100 [um], so pore diameters exceeding 100 um might overlap with neighbors.  The lower limit is to avoid vanishingly small pores.

>>> geom['pore.diameter'] = 0.00001 + sp.rand(pn.Np)*0.00099

This creates a ND-array of random numbers (between 0.00001 and 0.0001) that is *Np* long, meaning each pore is assigned a unique random number.

For throat diameter, we want them to always be smaller than the two pores which it connects to maintain physical consistency. This requires explaining how OpenPNM stores network topology.

>>> P12 = pn['throat.conns']  # An Nt x 2 list of pores on the end of each throat
>>> D12 = geom['pore.diameter'][P12]  # An Nt x 2 list of pore diameters
>>> Dt = sp.amin(D12, axis=1)  # An Nt x 1 list of the smaller pore from each pair
>>> geom['throat.diameter'] = Dt

Let's disect the above lines.  Firstly, P12 is a direct copy of the **Network's** 'throat.conns' array, which contains the indices of the pore pair connected by each throat.  Next, this *Nt-by-2* array is used to index into the 'pore.diameter' array, resulting in another *Nt-by-2* array containing the diameters of the pores connected by each throat.  Finally, the Scipy function ``amin`` is used to find the minimum diameter of each pore pair by specifying the ``axis`` keyword as 1, and the resulting *Nt-by-1* array is assigned to ``geom['throat.diameter']``.

Finally, we must specify the remaining geometrical properties of the pores and throats. Since we're creating a 'stick-and-ball' geometry, the sizes are calculated from the geometrical equations for spheres and cylinders.

For pore volumes, assume a sphere:

>>> Rp = geom['pore.diameter']/2
>>> geom['pore.volume'] = (4/3)*3.14159*(Rp)**3

The length of each throat is the center-to-center distance between pores, minus the radius of each of two neighbor pores.

>>> C2C = 0.0001  # The center-to-center distance between pores
>>> Rp12 = Rp[pn['throat.conns']]
>>> geom['throat.length'] = C2C - sp.sum(Rp12, axis=1)

The volume of each throat is found assuming a cylinder:

>>> Rt = geom['throat.diameter']/2
>>> Lt = geom['throat.length']
>>> geom['throat.volume'] = 3.14159*(Rt)**2*Lt

The basic geometrical properties of the network are now defined.

===============================================================================
Create Phases
===============================================================================

The simulation is now topologically and geometrically complete.  It has pore coordinates, pore and throat sizes and so on.  In order to perform any simulations it is necessary to define **Phase** objects that represent the fluids in the simulations:

>>> air = OpenPNM.Phases.GenericPhase(network=pn, name='air')
>>> water = OpenPNM.Phases.GenericPhase(network=pn, name='water')

``pn`` is passed as an argument because **Phases** must know to which **Network** they belong.  Also, note that ``pores`` and ``throats`` are NOT specified; this is because **Phases** are mobile and can exist anywhere or everywhere in the domain, so providing specific locations does not make sense.  Algorithms for dynamically determining actual phase distributions are discussed later.

.. note:: **Naming Objects**

	The above two lines also include a ``name`` argument.  All objects in OpenPNM can be named in this way if desired; however, if no name is given one will be generated.  The point of the name is to allow easy identification of an object at the command line, using the ``name`` attribute (``air.name``).  Objects can be renamed, so if you wish to override a default name simply use ``air.name = 'air'``.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Add Desired Properties to Phases
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Now it is necessary to fill these two **Phase** objects with the desired thermophysical properties.  The most basic means is to simply assign static values as follows:

>>> water['pore.temperature'] = 298.0
>>> water['pore.viscosity'] = 0.001
>>> air['pore.temperature'] = 298.0
>>> air['pore.viscosity'] = 0.0000173

OpenPNM includes a framework for calculating these type of properties from models and correlations, but this is beyond the aim of the present introductory tutorial.

.. note:: **Scalar to Vector Conversion During Assignment**

	The above lines illustrate a feature of OpenPNM that is worth pointing out now.  All pores need to have a diffusivity value associated with them; however, we often want to assign the same value to every pore.  If you assign a scalar value to any property in OpenPNM it will automatically be converted to a vector of the appropriate length (either *Np* or *Nt* long).  This is explained in more detail :ref:`here<inner_workings>`.

===============================================================================
Create Pore Scale Physics Objects
===============================================================================

We are still not ready to perform any simulations.  The last step is to define the desired pore scale physics models, which dictates how the phase and geometrical properties interact.  A classic example of this is the Hagen-Poiseuille equation for fluid flow through a throat, which predicts the flow rate as a function of the pressure drop  The flow rate is proportional to the geometrical size of the throat (radius and length) as well as properties of the fluid (viscosity).  It follows that this calculation needs to be performed once for each phase of interest since each has a different visocity.  This is accomlished by define a **Physics** object for each *Phase*:

>>> phys_water = OpenPNM.Physics.GenericPhysics(network=pn,
...                                             phase=water,
...                                             geometry=geom)
>>> phys_air = OpenPNM.Physics.GenericPhysics(network=pn,
...                                           phase=air,
...                                           geometry=geom)

**Physics** objects do not require the specification of which ``pores`` and ``throats`` where they apply, since this information is provided by the ``geometry`` argument which has already been assigned to specific locations.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Specify Desired Pore-Scale Models
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

We need to calculate the numerical values representing our chosen pore-scale physics.  To continue with the Hagen-Poiseuille example lets calculte the hydraulic conductance of each throat in the network.  The throat radius and length are easily accessed as:

>>> R = geom['throat.diameter']/2
>>> L = geom['throat.length']

The viscosity of the *Phases* was only defined in the pores; however, the hydraulic conductance must be calculated for each throat.  There are several options: (1) use a scalar value, (2) assign 'throat.viscosity' to each phase or (3) use interpolation to estimate throat viscosty as an average of the values in the neighboring pores.  The third option is suitable when there is a distribution of temperatures throughout the network and therefore visocity changes as well, and OpenPNM provides tools for this which are discussed later.  In the present case as simple scalar value is sufficient:

>>> mu_w = 0.001
>>> phys_water['throat.hydraulic_conductance'] = 3.14159*R**4/(8*mu_w*L)
>>> mu_a = 0.0000173
>>> phys_air['throat.hydraulic_conductance'] = 3.14159*R**4/(8*mu_a*L)

Note that both of these calcualation use the same geometrical properties (R and L) but different phase properties (mu_w and mu_a).

===============================================================================
Run Some Simulations
===============================================================================

Finally, it is now possible to run some simulations.  The code below estimates the permeabilty through the network by applying a pressure gradient across and calculating the flux.  This starts by creating a StokesFlow *Algorithm*, which is pre-defined in OpenPNM:

>>> alg = OpenPNM.Algorithms.StokesFlow(network=pn, phase=air)

Like all the above objects, algorithms must be assigned to a *Network* via the ``network`` argument.  This algorithm is also associated with a *Phase* object, in this case ``air``, which dictates which pore-scale *Physics* properties to use (recall that ``phys_air`` was associated with ``air``).

Next the boundary conditions are applied using the ``set_boundary_conditions`` method on the *Algorithm* object.  Let's apply a 1 atm pressure gradient between the left and right sides of the domain:

>>> BC1_pores = pn.pores('right')
>>> alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=202650, pores=BC1_pores)
>>> BC2_pores = pn.pores('left')
>>> alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=101325, pores=BC2_pores)

To actually run the algorithm use the ``run`` method.  This builds the coefficient matrix from the existing values of hydraulic conductance, and inverts the matrix to solve for pressure in each pore, and stores the results within the *Algorithm's* dictionary under 'pore.pressure'\:

>>> alg.run()

The results ('pore.pressure') are held within the ``alg`` object and must be explicitly returned to the ``air`` object by the user if they wish to use these values in a subsequent calcualation.  The point of this data containment is to prevent unwanted overwriting of data.  Each algorithm has a method called ``return_results`` which places the pertinent values back onto the appropriate *Phase* object.

>>> alg.return_results()

===============================================================================
Visualise the Results
===============================================================================
We can now visualise our network and simulation results.  OpenPNM does not support native visualization, so data must be exported to a file for exploration in another program such as any of the several VTK front ends (i.e. Paraview).

>>> OpenPNM.export(network=pn, filename='net.vtp')

This creates a *net.vtp* file in the active directory, which can be loaded from ParaView. For a quick tutorial on the use of Paraview with OpenPNM data, see :ref:`Using Paraview<paraview_example>`.
