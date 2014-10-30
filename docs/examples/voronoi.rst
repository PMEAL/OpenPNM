.. _gostick:

###############################################################################
Random PNM with Delaunay Network and Voronoi Geometry
###############################################################################

.. _J.T. Gostick et al. / JES 160 (2013) F731-F743: http://jes.ecsdl.org/cgi/doi/10.1149/2.009308jes

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Getting Started
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

In this tutorial, we will demonstrate concepts of random pore network modelling outlined in J.T. Gostick's 2013 paper `[1]`_. Pores are randomly dispersed inside the domain and connections are defined by the nearest neighbour Delaunay triangulation algorithm. The Voronoi diagram is the copliment to the Delaunay trinagulation and is used to replicate the fibres in a porous domain such as a fuel cell GDL. Each pore is surrounded by a cage of fibres which are created from the intersections of the equidistant planes between neighbouring pores. Where the planes intersect, a Voronoi vertex is created, and these are saved to the Delaunay Network object as "pore.vertices". A throat is defined by a set of vertices which are shared by two neighbouring pores. The throat vertices will be coplanar and the throat normal is the vector from one pore's co-ordinate to the other. The vertices are used by the Voronoi geometry model which effectively erodes the fibre space by offsetting the vertices in the planes of each connecting throat by a fibre radius parameter passed to the method.
 
.. _[1]: http://jes.ecsdl.org/cgi/doi/10.1149/2.009308jes

There will be a general layout to complete this simulation: 

1. Set up network 
2. Set up geometry and geometrical methods and inspect the pore and throat geometry
3. Set up phases and methods 
4. Set up phase physics and methods 
5. Run ordinary percolation 
6. Plot a drainage curve

We first import the OpenPNM code.

.. code-block:: python
    
    import OpenPNM
   
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Setting up Network and Geometry
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

To begin our simulation, we must first generate our random network and geometry.  This includes:

1. Creating a Delaunay network object and a Voronoi geometry object
2. Creating boundary pores
3. Assingning a 'Voronoi' geometry to the initial set of pores and a 'Boundary' geometry to the boundary pores

.. code-block:: python

	pn = OpenPNM.Network.Delaunay(num_pores=100, domain_size=[1,1,1],name='net')
	pn.add_boundaries()
	# create geometry
	fibre_rad=0.03
	Ps = pn.pores('boundary',mode='not')
	Ts = pn.find_neighbor_throats(pores=Ps,mode='intersection',flatten=True)
	geo = OpenPNM.Geometry.Voronoi(network=pn,pores=Ps,throats=Ts,fibre_rad=fibre_rad)
	# create boundary geometry
	Ps = pn.pores('boundary')
	Ts = pn.find_neighbor_throats(pores=Ps,mode='not_intersection')
	boun = OpenPNM.Geometry.Boundary(network=pn,pores=Ps,throats=Ts,name='boun')
	
It is worth mentioning a little about the boundaries at this point. Behind the scenes all the pores in the network were temporarily reflected about the planes confining the pore space during the Delaunay tesselation. This has the effect of creating throats on the outer confines that are aligned with the boundary planes. The add_boundaries() method makes use of this feature and generates 6 new sets of boundary pores which have their co-ordinates placed on the boundary plane in-line with the adjoining pore. The boundary pores are labelled "left_boundary", "right_boundary", etc. and all boundary pores are also given the label "boundary". It is generally a good idea to assign these pores to a Boundary type Geometry which will give the pores zero volume.
	
The 'Voronoi' geometry object has two visualisation methods which are unique to it. Both methods will plot a selection of pores or throats and help visualise and check the offsetting routine.  

.. code-block:: python

	#print all pores in the Voronoi geometry
	geo.print_pore(geo['pore.map'])
	throats = pn.find_neighbor_throats(pores=[0])
	#print all throats connected to the first pore in the network
	geo.print_throat(throats)

.. image:: http://imgur.com/zrEFZ1O.png

.. image:: http://imgur.com/mmnQsHz.png

We can also see the distributions of the pore and throat sizes.

.. code-block:: python

	OpenPNM.Postprocessing.Plots.distributions(net=pn)

.. image:: http://imgur.com/Zw93uSV.png

Three new pieces of information are added to the Network object when the Vornoi geometry is applied to a set of pores: 'throat.verts', 'throat.offset_verts' and 'throat.normals'. The offset vertices are subsequently used to calculate the throat area. Also the pore volume is calculated using the convex hull of all the offset throat vertices for a given pore.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Setting up the Phases and Physics
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Now we are ready to set up our phases (water and air) and the physics corresponding to each of these phases. OpenPNM has built in air and water phases, so we can use those.

.. code-block:: python

    #set up phases
    air = OpenPNM.Phases.Air(network = pn, name = 'air')
    water = OpenPNM.Phases.Water(network = pn, name = 'water')

We are now ready to establish physical properties for our fluid objects. To do this, we will: 

1. Create physics objects associated with our fluids (by using BasePhyics we don't have to add methods for calculating each property because they are already included) 
2. Use our regenerate_physics() method to calculate these properties

.. code-block:: python

    #create physics objects associated with our phases
    Ps = pn.pores()
    Ts = pn.throats()
    phys_water = OpenPNM.Physics.Standard(network=pn,phase=water,pores=Ps,throats=Ts,dynamic_data=True,name='standard_water_physics')
    phys_air = OpenPNM.Physics.Standard(network=pn,phase=air,pores=Ps,throats=Ts,dynamic_data=True,name='standard_air_physics')
	
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Running Ordinary Percolation & Visualising the Output
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

A simple algorithm to demonstrate the features of the network is the Ordinary Percolation algorithm.  
To run our simulation, we will:

1. Pick inlet pores
2. Create an Ordinary Percolation algorithm object
3. Run our algorithm object
4. Call update() so that occupancy of pores and throats for each fluid will be set

.. code-block:: python

    inlets = pn.pores('bottom_boundary')
    used_inlets = [inlets[x] for x in range(0, len(inlets), 2)]
    
    #using every other pore in the bottom and boundary as an inlet
    #prevents extremely small diffusivity and permeability values in the z direction
    used_inlets = [inlets[x] for x in range(0, len(inlets), 2)]
    
    OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=pn,invading_phase=water,defending_phase=air)
    OP_1.run(inlets=used_inlets,npts=100)

This algorithm performed a start to finish simulation, which fully flooded the network. The 'return_results()' command can be used to update the phase occupancy values throughout the network. 

.. code-block:: python

	#Update the simulation until saturation is at 50%
	OP_1.return_results(sat=0.5)

OpenPNM makes it very easy to inspect the ouput of the algorithm through the "Postprocessing" methods.  

.. code-block:: python

	OpenPNM.Postprocessing.Plots.drainage_curves(OP_1,timing=None)

We can also view the network data by creating vtk files to be viewed using ParaView (downloadable at http://www.paraview.org/download/ ). It is suggested that version 3.98 is downloaded instead of 4.1).  If we visualize our pore network model with phase data included it will look like this:

.. image:: http://imgur.com/lmjSHG7.png
	
Spherical glyphs are used to represent the pores and are sized using the pore diameter. The water.occupancy data is used to colour the glyphs and those that are un-occupied are set to be invisible using the opacity scale.

To create the vtk file use the following command

.. code-block:: python
	
    import OpenPNM.Utilities.IO as io
    io.VTK.save(network=pn,phases=[air,water])

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Differences between OpenPNM and Gostick's simulation
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

The major difference between the technique described in Gostick's paper `[1]`_ and that implemented in OpenPNM is the method used to replicate the fibres. Where Gostick uses morphological image analysis to construct 3D voxel images of the solid fibres a simpler method using geometrical offsetting of the Voronoi vertices is used in OpenPNM. This technique is quicker than image analysis and can be achieved entirely within the framework of OpenPNM, however, a fibrous pore space is less well represented as fibres are effectively treated as polygons rather than cylinders.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
References
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

`[1]`_ J. T. Gostick, "Random Pore Network Modeling of Fibrous PEMFC Gas Diffusion Media Using Voronoi and Delaunay Tessellations" Journal of the Electrochemical Society, vol. 160, issue 8, pp. F731-F743, 2013.