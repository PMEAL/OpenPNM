.. _gostick:

###############################################################################
Random PNM with Delaunay Network and Voronoi Geometry
###############################################################################

.. _J.T. Gostick et al. / JES 160 (2013) F731-F743: http://jes.ecsdl.org/cgi/doi/10.1149/2.009308jes

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Getting Started
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

In this tutorial, we will demonstrate concepts of random pore network modelling outlined in J.T. Gostick's 2013 paper `[1]`_. Pores are randomly dispersed inside the domain and connections are defined by the nearest neighbour Delaunay triangulation algorithm. The Voronoi diagram is the compliment to the Delaunay trinagulation and is used to replicate the fibres in a porous domain such as a fuel cell GDL. Each pore is surrounded by a cage of fibres which are created from the intersections of the equidistant planes between neighbouring pores. Where the planes intersect, a Voronoi vertex is created, and these are saved to the Delaunay Network object as "pore.vert_index" and "throat.vert_index". A throat is defined by a set of vertices which are shared by two neighbouring pores. The throat vertices will be coplanar and the throat normal is the vector from one pore's co-ordinate to the other. The vertices are used by the Voronoi geometry model which creates a 3D image of the fibres using a supplied fibre radius. Image analysis is then performed to extract pore and throat sizes by using the convex hull of the pore and throat vertices.

.. _[1]: http://jes.ecsdl.org/cgi/doi/10.1149/2.009308jes

There will be a general layout to complete this simulation:

1. Set up network
2. Set up geometry and geometrical methods and inspect the pore and throat geometry
3. Set up phases and methods
4. Set up phase physics and methods
5. Run ordinary percolation
6. Plot a drainage curve

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Setting up Network and Geometry
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

We first import the OpenPNM code.

.. code-block:: python

    import OpenPNM

To begin our simulation, we must first generate our random network and geometry.  This includes:

1. Getting a handle to the OpenPNM Controller and setting the loglevel to display useful information
2. Creating a Delaunay network object and a Voronoi geometry object
3. Assingning a 'Voronoi' geometry to the initial set of pores and a 'Boundary' geometry to the boundary pores

.. code-block:: python

	ctrl = OpenPNM.Base.Controller()
	ctrl.loglevel = 20
	pn = OpenPNM.Network.Delaunay(num_pores=100, domain_size=[1e-4, 1e-4, 1e-4], name='net')
	pn.add_boundaries()
	# create geometry
	fibre_rad = 5e-6
	Ps = pn.pores()
	Ts = pn.find_neighbor_throats(pores=Ps, mode='intersection', flatten=True)
	geo = OpenPNM.Geometry.Voronoi(network=pn, pores=Ps, throats=Ts, fibre_rad=fibre_rad, voxel_vol=True, vox_len=1e-6, name='vor')

It is worth mentioning a little about the boundaries at this point. Behind the scenes all the pores in the network were temporarily reflected about the planes confining the pore space during the Delaunay tesselation. This has the effect of creating throats on the outer confines that are aligned with the boundary planes. The add_boundaries() method makes use of this feature and generates 6 new sets of boundary pores which have their co-ordinates placed on the boundary plane in-line with the adjoining pore. The boundary pores are labelled "left_boundary", "right_boundary", etc. and all boundary pores are also given the label "boundary".

The 'vertexops' utility has two visualisation methods which are unique to Voronoi geometries. Both methods will plot a selection of pores or throats and help visualise and check the offsetting routine.

.. code-block:: python

	#plot all pores in the Voronoi geometry
	import OpenPNM.Utilities.vertexops as vo
	vo.plot_pore(geo, geo.pores())
	throats = pn.find_neighbor_throats(pores=[0])
	#plot all throats connected to the first pore in the network
	vo.plot_throat(geo, throats)

.. image:: http://imgur.com/icZfxQ2.png

.. image:: http://imgur.com/sCL1kY1.png

We can also see the distributions of the pore and throat sizes.

.. code-block:: python

	OpenPNM.Postprocessing.Plots.distributions(pn)

.. image:: http://imgur.com/IylwVk5.png

The Voronoi geometry has two parameters which control whether a voxel image of the fibres is generated or not and what resolution is used. If 'voxel_vol' is set to True then the image is created with resolution set by the 'vox_len' parameter. The default is to generate the image at resolution of 1e-6 voxel side length giving each voxel a volume of 1e-18. Generating the voxel image is a memory intensive process relying on many image analysis routines and it is recommended that a smaller network is tested first on your machine whilst monitroing your system performance to guage whether larger networks are possible. Setting 'voxel_vol' to False will mean that pore volumes are calculated from the ConvexHull of the offset throat vertices which is a faster but less accurate method.

If the voxel image is generated a few more methods of the Voronoi Geometry class can be used to visualize and export the image.

.. code-block:: python

	geo.plot_fibre_slice(plane=[0.5,0,0])

.. image:: http://imgur.com/GMzeBck.png

.. code-block:: python

	geo.plot_porosity_profile()

.. image:: http://imgur.com/j0kTCRj.png

The following command creates a binary array which can be used in Matlab.

.. code-block:: python

	geo._export_fibre_image()

In addition to visualization within OpenPNM using the voxel image it is possible to export a pickle dump of the throat vertices in convex hull order which form the skeleton of the fibres.

.. code-block:: python

	pn._export_vor_fibres()

This can then be used in Blender to create images such as this:

.. image:: http://imgur.com/MCPQ1S7.png

The Blender script used to create the multi-faceted Tetrahedron is as follows:

.. code-block:: python

	import bpy
	import mathutils
	import numpy as np
	import pickle as pickle

	class TetrahedronMakerPanel(bpy.types.Panel):
		bl_space_type = "VIEW_3D"
		bl_region_type = "TOOLS"
		bl_context = "objectmode"
		bl_category = "Create"
		bl_label = "Add Tetrahedron"

		def draw(self, context):
			TheCol = self.layout.column(align=True)
			TheCol.operator("mesh.make_tetrahedron", text="Add Tetrahedron")
		#end draw

	#end TetrahedronMakerPanel

	class MakeTetrahedron(bpy.types.Operator):
		bl_idname = "mesh.make_tetrahedron"
		bl_label = "Add Tetrahedron"
		bl_options = {"UNDO"}

		def invoke(self, context, event):
			data = pickle.load( open( "C:\\Working_Files\\OpenPNM\\fibres.p", "rb" ) )
			Verts = np.around(data["Verts"]*1e4,2)
			Faces = data['Indices']
			NewMesh = bpy.data.meshes.new("Tetrahedron")
			NewMesh.from_pydata \
			  (
				Verts,
				[],
				Faces
			  )
			NewMesh.update()
			NewObj = bpy.data.objects.new("Tetrahedron", NewMesh)
			context.scene.objects.link(NewObj)
			return {"FINISHED"}
		#end invoke

	#end MakeTetrahedron

	bpy.utils.register_class(MakeTetrahedron)
	bpy.utils.register_class(TetrahedronMakerPanel)

The image is then created by deleteing the faces leaving just the wire-frame edges of the Voronoi diagram and then converting the object to a curve and applying a circular bevel curve.

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Setting up the Phases and Physics
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Now we are ready to set up our phases (water and air) and the physics corresponding to each of these phases. OpenPNM has built in air and water phases, so we can use those.

.. code-block:: python

    #set up phases
    air = OpenPNM.Phases.Air(network = pn, name = 'air')
    water = OpenPNM.Phases.Water(network = pn, name = 'water')

We are now ready to establish physical properties for our fluid objects. To do this, we will create physics objects associated with our fluids (by using OpenPNM.Physics.Standard we don't have to add methods for calculating each property because they are already included)

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
4. Return results so that occupancy of pores and throats for each fluid will be set

.. code-block:: python

    inlets = pn.pores('bottom_boundary')
    used_inlets = [inlets[x] for x in range(0, len(inlets), 2)]

    #using every other pore in the bottom and boundary as an inlet
    #prevents extremely small diffusivity and permeability values in the z direction
    used_inlets = [inlets[x] for x in range(0, len(inlets), 2)]

    OP_1 = OpenPNM.Algorithms.OrdinaryPercolation(network=pn, invading_phase=water, defending_phase=air)
    OP_1.run(inlets=used_inlets, npts=100)

This algorithm performed a start to finish simulation, which fully flooded the network. The 'return_results()' command can be used to update the phase occupancy values throughout the network.

.. code-block:: python

	#Update the simulation until saturation is at 50%
	OP_1.return_results(sat=0.5)

OpenPNM makes it very easy to inspect the ouput of the algorithm through the "Postprocessing" methods.

.. code-block:: python

    import matplotlib.pyplot as plt
    fig=plt.figure()
    OpenPNM.Postprocessing.Plots.drainage_curves(OP_1, timing=None, fig=fig)

We can also view the network data by creating vtk files to be viewed using ParaView (downloadable at http://www.paraview.org/download/ ). It is suggested that version 3.98 is downloaded instead of 4.1).  If we visualize our pore network model with phase data included it will look like this:

.. image:: http://imgur.com/lmjSHG7.png

Spherical glyphs are used to represent the pores and are sized using the pore diameter. The water.occupancy data is used to colour the glyphs and those that are un-occupied are set to be invisible using the opacity scale.

To create the vtk file use the following command

.. code-block:: python

    ctrl.export(pn)

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
References
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

`[1]`_ J. T. Gostick, "Random Pore Network Modeling of Fibrous PEMFC Gas Diffusion Media Using Voronoi and Delaunay Tessellations" Journal of the Electrochemical Society, vol. 160, issue 8, pp. F731-F743, 2013.
