.. _wu:

###############################################################################
Example: Regenerating Data from `R. Wu et al. / Elec Acta v54n25 (2010) 7394–7403`_
###############################################################################

.. _R. Wu et al. / Elec Acta v54n25 (2010) 7394–7403: http://www.sciencedirect.com/science/article/pii/S0013468610009503

.. warning::

    The following example was painstakingly created using an early version OpenPNM, and does not work properly with the current version.  It still valuable as it provides a detailed illustration of the OpenPNM workflow.  

.. code-block:: python

	import OpenPNM
	import matplotlib.pyplot as plt
	import scipy as sp
	import numpy as np

	n = 8
	Lc = 25e-6

	pn = OpenPNM.Network.Cubic(name = 'Wu')

	#code to run if we want to set add_boundaries to be False
	pn.generate(divisions = [n,n,2*n], add_boundaries= False, lattice_spacing = [Lc])

	#code to run if boundaries was set to false
	Ps = pn.pores()
	Ts = pn.throats()
	geo = OpenPNM.Geometry.Toray090(name='wu_geometry',network=pn,pores=Ps,throats=Ts)

	#code to run if boundaries was set to True
	#pn.generate(divisions = [n,n,2*n], add_boundaries= True, lattice_spacing = [Lc])
	#geo = OpenPNM.Geometry.GenericGeometry(name = 'wu_geometry', network = pn)
	#geo.set_locations(pores = pn.pores('internal'), throats = 'all')
	#boun = pn.add_geometry(subclass='Boundary',name='boun')
	#boun.set_locations(pores=pn.pores('boundary'))

	low = .5e-6
	high = 9.5e-6

	#geo.add_prop(prop='pore_diameter',model='constant', value = 24e-6)
	#geo.add_method(prop='pore_volume',model='sphere')
	#geo.add_method(prop='throat_length',model='straight')
	#geo.add_method(prop='throat_volume',model='cylinder')
	geo['pore.diameter'] = 24e-6

	#setting throat diameters to be a uniform distribution
	radii = low + sp.rand(pn.num_throats())*(high - low)
	geo.set_data(prop = 'diameter', throats = pn.throats(), data = radii*2)

	#geo.regenerate()


	air = OpenPNM.Phases.Air(network = pn, name = 'air')
	water = OpenPNM.Phases.Water(network = pn, name = 'water')
	air.regenerate()
	water.regenerate()


	def bulk_diffusion_wu(physics,
						  network,
						  phase,
						  geometry,
						  propname,
						  diffusivity = 'diffusivity',
						  molar_density = 'molar_density',
						  throat_diameter = 'diameter',
						  throat_length = 'length',
						  pore_diameter = 'diameter',
						  **params):
		r"""
			Calculate the diffusive conductance of throats in network (instead of a
			conduit) based on the areas

			Parameters
			----------
			network : OpenPNM Network Object

			phase : OpenPNM Phase Object
			The phase of interest

			Notes
			-----
			This function requires that all the necessary phase properties already be
			calculated.

			"""
		#ct = phase.get_data(prop='molar_density',throats='all',mode='interpolate')
		#Interpolate pore values to throats
		DABt = phase.get_data(prop='diffusivity',throats=geometry.throats(),mode='interpolate')
		#Find g for full throat
		tdia = network.get_throat_data(prop=throat_diameter)
		tlen = network.get_throat_data(prop=throat_length)
		gt = (sp.pi*DABt*tdia**2)/(tlen*4)
		g = gt[geometry.throats()]
		phase.set_data(prop=propname,throats=geometry.throats(),data=g)

	Ps = geo.pores()
	Ts = geo.throats()
	phys_water = OpenPNM.Physics.Standard(network=pn,phase=water, pores=Ps, throats=Ts, name='standard_water_physics')
	phys_air = OpenPNM.Physics.Standard(network=pn,phase=air, pores=Ps, throats=Ts, geometry = geo, name='standard_air_physics')

	#phys_water.add_model(prop='capillary_pressure', model='washburn') #accounts for cylindrical throats
	#phys_water.add_model(prop='hydraulic_conductance',model='hagen_poiseuille')
	#phys_water.add_model(prop='diffusive_conductance', model='bulk_diffusion', shape = 'circular')
	#phys_air.add_model(prop='hydraulic_conductance',model='hagen_poiseuille')

	bulk_diffusion_wu(physics = phys_air, network = pn, phase = air, geometry = geo, propname = 'diffusive_conductance')
	phys_water.regenerate()
	phys_air.regenerate()

	inlets = pn.get_pore_indices(labels = ['bottom']) #can put in brackets so the whole bottom of the lattice is considered 1 inlet
	outlets = pn.get_pore_indices(labels = ['top'])

	IP_1 = OpenPNM.Algorithms.InvasionPercolation(network = pn, name = 'OP_1')
	IP_1.setup(invading_phase = water, defending_phase = air, inlets = inlets, outlets = outlets, end_condition = 'total')
	IP_1.run()

	max_inv_seq = max(IP_1.get_pore_data(prop = 'IP_inv_seq'))
	x_values = []
	y_values = []

	for x in range(50):
		IP_1.return_results(IPseq = max_inv_seq*(x/50.0))

		phys_air.add_model(model=OpenPNM.Physics.models.multiphase.conduit_conductance,
				   propname='throat.conduit_diffusive_conductance',
				   throat_conductance='throat.diffusive_conductance',
				   mode='strict')
		phys_water.add_model(model=OpenPNM.Physics.models.multiphase.conduit_conductance,
				   propname='throat.conduit_diffusive_conductance',
				   throat_conductance='throat.diffusive_conductance',
				   mode='strict')
		phys_air.add_model(model=OpenPNM.Physics.models.multiphase.conduit_conductance,
				   propname='throat.conduit_hydraulic_conductance',
				   throat_conductance='throat.hydraulic_conductance',
				   mode='strict')
		phys_water.add_model(model=OpenPNM.Physics.models.multiphase.conduit_conductance,
				   propname='throat.conduit_hydraulic_conductance',
				   throat_conductance='throat.hydraulic_conductance',
				   mode='strict')

		Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(loggername = 'Fickian', name = 'Fickian', network = pn)

		#set labels for top boundary
		#set labels for bottom boundary
		A = pn._Nx**2
		top_pores = pn.get_pore_indices('top')
		bottom_pores = pn.get_pore_indices('bottom')
		z_dimension = int(pn.domain_length(top_pores,bottom_pores)/Lc) #number of pores in the z direction
		quarter_layer = z_dimension/4 #estimates which layer marks 1/4 up the lattice
		pore_number = int(quarter_layer*A) #gives the first pore in the layer 1/4 up the lattice

		bottom_boundary = list(range(pore_number, pore_number + A))
		top_boundary = list(range(pn.num_pores() - pore_number, pn.num_pores() - pore_number +A))

		Fickian_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=bottom_boundary)
		Fickian_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.2, pores=top_boundary)

		Fickian_alg.setup(conductance = 'conduit_diffusive_conductance',phase=air)
		Fickian_alg.run()

		effective_diffusivity = Fickian_alg.calc_eff_diffusivity(clean = False)
		bulk_diffusivity = air.get_pore_data(prop = 'diffusivity')

		#calculation of saturation
		center_pores = list(range(bottom_boundary[-1], top_boundary[0]))
		final_pores = water['pore.occupancy'][center_pores]
		pore_volumes = pn['pore.volume'][center_pores]

		saturation = sum(final_pores*pore_volumes)/sum(pore_volumes)

		x_values.append(saturation)
		y_values.append((effective_diffusivity/bulk_diffusivity)[0])


	plt.plot(x_values, y_values, 'ro')
	plt.title('normalized diffusivity versus saturation')
	plt.xlabel('saturation')
	plt.ylabel('normalized diffusivity')
	plt.show()

	x_values = []
	y_values = []

	for x in range(20):
		n = 8
		Lc = 25e-6
		np.random.seed()

		pn = OpenPNM.Network.Cubic(name = 'Wu')

		#code to run if we want to set add_boundaries to be False
		pn.generate(divisions = [n,n,2*n], add_boundaries= False, lattice_spacing = [Lc])

		#code to run if boundaries was set to false
		Ps = pn.pores()
		Ts = pn.throats()
		geo = OpenPNM.Geometry.Toray090(name='wu_geometry',network=pn,pores=Ps,throats=Ts)

		low = .5e-6
		high = 9.5e-6

		geo['pore.diameter'] = 24e-6

		#setting throat diameters to be a uniform distribution
		radii = low + np.random.random(pn.num_throats())*(high - low)
		geo.set_data(prop = 'diameter', throats = pn.throats(), data = radii*2)

	#    geo.regenerate()

		air = OpenPNM.Phases.Air(network = pn, name = 'air')
		water = OpenPNM.Phases.Water(network = pn, name = 'water')
		air.regenerate()
		water.regenerate()

		Ps = geo.pores()
		Ts = geo.throats()
		phys_water = OpenPNM.Physics.Standard(network=pn,phase=water, pores=Ps, throats=Ts, name='standard_water_physics')
		phys_air = OpenPNM.Physics.Standard(network=pn,phase=air, pores=Ps, throats=Ts, geometry = geo, name='standard_air_physics')

	#    phys_water.add_method(prop='capillary_pressure', model='washburn') #accounts for cylindrical throats
	#    phys_water.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
	#    phys_water.add_method(prop='diffusive_conductance', model='bulk_diffusion', shape = 'circular')
	#    phys_air.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')

		bulk_diffusion_wu(physics = phys_air, network = pn, phase = air, geometry = geo, propname = 'diffusive_conductance')
		phys_water.regenerate()
		phys_air.regenerate()

		inlets = pn.get_pore_indices(labels = ['bottom']) #can put in brackets so the whole bottom of the lattice is considered 1 inlet
		outlets = pn.get_pore_indices(labels = ['top'])

		IP_1 = OpenPNM.Algorithms.InvasionPercolation(network = pn, name = 'OP_1')
		IP_1.setup(invading_phase = water, defending_phase = air, inlets = inlets, outlets = outlets, end_condition = 'total')
		IP_1.run()

		max_inv_seq = max(IP_1.get_pore_data(prop = 'IP_inv_seq'))

		for x in range(50):
			IP_1.return_results(IPseq = max_inv_seq*(x/50.0))

			phys_air.add_model(model=OpenPNM.Physics.models.multiphase.conduit_conductance,
					   propname='throat.conduit_diffusive_conductance',
					   throat_conductance='throat.diffusive_conductance',
					   mode='strict')
			phys_water.add_model(model=OpenPNM.Physics.models.multiphase.conduit_conductance,
					   propname='throat.conduit_diffusive_conductance',
					   throat_conductance='throat.diffusive_conductance',
					   mode='strict')
			phys_air.add_model(model=OpenPNM.Physics.models.multiphase.conduit_conductance,
					   propname='throat.conduit_hydraulic_conductance',
					   throat_conductance='throat.hydraulic_conductance',
					   mode='strict')
			phys_water.add_model(model=OpenPNM.Physics.models.multiphase.conduit_conductance,
					   propname='throat.conduit_hydraulic_conductance',
					   throat_conductance='throat.hydraulic_conductance',
					   mode='strict')

			Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(loggername = 'Fickian', name = 'Fickian', network = pn)

			#set labels for top boundary
			#set labels for bottom boundary
			A = pn._Nx**2

			top_pores = pn.get_pore_indices('top')
			bottom_pores = pn.get_pore_indices('bottom')
			z_dimension = int(pn.domain_length(top_pores,bottom_pores)/Lc) #number of pores in the z direction
			quarter_layer = z_dimension/4 #estimates which layer marks 1/4 up the lattice
			pore_number = int(quarter_layer*A) #gives the first pore in the layer 1/4 up the lattice

			bottom_boundary = list(range(pore_number, pore_number + A))
			top_boundary = list(range(pn.num_pores() - pore_number, pn.num_pores() - pore_number +A))

			Fickian_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=bottom_boundary)
			Fickian_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.2, pores=top_boundary)

			Fickian_alg.setup(conductance = 'conduit_diffusive_conductance',phase=air)
			Fickian_alg.run()

			effective_diffusivity = Fickian_alg.calc_eff_diffusivity(clean = False)
			bulk_diffusivity = air.get_pore_data(prop = 'diffusivity')

			#calculation of saturation
			center_pores = list(range(bottom_boundary[-1], top_boundary[0]))
			final_pores = water['pore.occupancy'][center_pores]
			pore_volumes = pn['pore.volume'][center_pores]

			saturation = sum(final_pores*pore_volumes)/sum(pore_volumes)

			x_values.append(saturation)
			y_values.append((effective_diffusivity/bulk_diffusivity)[0])

	from matplotlib.font_manager import FontProperties
	fontP = FontProperties()
	fontP.set_size('small')

	wu_average_x_values = [0.004, 0.021, 0.052, 0.081, 0.129, 0.162, 0.186, 0.219, 0.261,
						   0.286, 0.324, 0.363, 0.42, 0.478, 0.531, 0.586, 0.64, 0.698, 0.747, 0.802]
	wu_average_y_values = [0.118, 0.113, 0.105, 0.096, 0.085, 0.078, 0.07, 0.062, 0.054, 0.049, 0.04,
						   0.033, 0.027, 0.02, 0.012, 0.006, 0.003, 0.002, 0.002, 0.002]

	p1, = plt.plot(x_values, y_values, 'wo')
	p2, = plt.plot(wu_average_x_values, wu_average_y_values, 'ro')
	plt.title('normalized diffusivity versus saturation')
	plt.xlabel('saturation')
	plt.ylabel(r'$\frac{D_e}{D_b}$')
	plt.ylim([0, .15])
	plt.xlim([0, 1])
	plt.legend([p1, p2],
			   [r'$\frac{D_e}{D_b} = f(\epsilon, \phi)g(s, \phi)$' + '\n' + r'$X = 1.8$' +
			   '\n' + r'$Z_t = 2.0$' + '\n' + r'$Z_i = 4.0$' + '\n' + r'$\beta = 1.0$' + '\n' + r'$n = 14$', "Wu's results"])

	plt.show()

	y_2_values = []

	n_values = [8, 10, 12, 14, 16, 18, 20]

	for x in range(5):
		for n in n_values:

			Lc = 25e-6

			pn = OpenPNM.Network.Cubic(name = 'Wu')

			#code to run if we want to set add_boundaries to be False
			pn.generate(divisions = [n,n,2*n], add_boundaries= False, lattice_spacing = [Lc])

			Ps = pn.pores()
			Ts = pn.throats()
			geo = OpenPNM.Geometry.Toray090(name='wu_geometry',network=pn,pores=Ps,throats=Ts)

			low = .5e-6
			high = 9.5e-6

			geo.set_data(prop = 'diameter', pores = Ps, data = 24e-6, mode = 'overwrite')

			#setting throat diameters to be a uniform distribution
			radii = low + sp.rand(pn.num_throats())*(high - low)
			geo.set_data(prop = 'diameter', throats = pn.throats(), data = radii*2)
	#
	#        geo.regenerate()

			#phases
			air = OpenPNM.Phases.Air(network = pn, name = 'air')
			water = OpenPNM.Phases.Water(network = pn, name = 'water')
			air.regenerate()
			water.regenerate()

			#physics objects
			phys_water = OpenPNM.Physics.GenericPhysics(network=pn,phase=water, geometry = geo, name='standard_water_physics')
			phys_air = OpenPNM.Physics.GenericPhysics(network=pn,phase=air, geometry = geo, name='standard_air_physics')

	#        phys_water.add_method(prop='capillary_pressure', model='washburn') #accounts for cylindrical throats
	#        phys_water.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')
	#        phys_water.add_method(prop='diffusive_conductance', model='bulk_diffusion', shape = 'circular')
	#        phys_air.add_method(prop='hydraulic_conductance',model='hagen_poiseuille')

			bulk_diffusion_wu(physics = phys_air, network = pn, phase = air, geometry = geo, propname = 'diffusive_conductance')
			phys_water.regenerate()
			phys_air.regenerate()

			#Invasion percolation
			inlets = pn.get_pore_indices(labels = ['bottom']) #can put in brackets so the whole bottom of the lattice is considered 1 inlet
			outlets = pn.get_pore_indices(labels = ['top'])

			air.set_data(pores = pn.pores(), prop = 'occupancy', data = 1)
			air.set_data(throats = pn.throats(), prop = 'occupancy', data = 1)

			water.set_data(pores = pn.pores(), prop = 'occupancy', data = 0)
			water.set_data(throats = pn.throats(), prop = 'occupancy', data = 0)

			Fickian_alg = OpenPNM.Algorithms.FickianDiffusion(loggername = 'Fickian', name = 'Fickian', network = pn)

			#set labels for top boundary
			#set labels for bottom boundary
			A = pn._Nx**2

			top_pores = pn.get_pore_indices('top')
			bottom_pores = pn.get_pore_indices('bottom')
			z_dimension = int(pn.domain_length(top_pores,bottom_pores)/Lc) #number of pores in the z direction
			quarter_layer = z_dimension/4 #estimates which layer marks 1/4 up the lattice
			pore_number = int(quarter_layer*A) #gives the first pore in the layer 1/4 up the lattice

			bottom_boundary = list(range(pore_number, pore_number + A))
			top_boundary = list(range(pn.num_pores() - pore_number, pn.num_pores() - pore_number +A))

			Fickian_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.6, pores=bottom_boundary)
			Fickian_alg.set_boundary_conditions(bctype='Dirichlet', bcvalue=0.2, pores=top_boundary)

			Fickian_alg.setup(conductance = 'diffusive_conductance',phase=air)
			Fickian_alg.run()

			effective_diffusivity = Fickian_alg.calc_eff_diffusivity(clean = False)
			bulk_diffusivity = air.get_pore_data(prop = 'diffusivity')

			y_2_values.append((effective_diffusivity/bulk_diffusivity)[0])



	plt.plot(n_values + n_values + n_values + n_values + n_values, y_2_values, 'ro')
	plt.title('F(epsilon) versus N')
	plt.xlabel('N')
	plt.ylabel('F(epsilon)')
	plt.axis(xmin = 6,xmax = 22,ymin= 0,ymax = .2)
	plt.show()


	normalize_factor = y_values[0]
	g_values = list(range(len(y_values)))
	for x in range(1000):
		if x%50 == 0:
			normalize_factor = y_values[x]
		g_values[x] = y_values[x] / normalize_factor

	wu_saturation = [0.004, 0.066, 0.0930, .119, 0.14, 0.175, 0.209, 0.24, 0.282, 0.32, 0.371, 0.413,
		0.464, 0.517, 0.605, 0.672, 0.761, 0.831, 0.898, 0.948, 0.996]
	wu_g_values = [0.986, 0.838, 0.758, 0.701, 0.651, 0.576, 0.516, 0.456, 0.39, 0.335, 0.268, 0.221,
		0.171, 0.111, 0.067, 0.04, 0.019, 0.007, 0.003, 0.003, 0.003]

	p1, = plt.plot(x_values, g_values, 'wo')
	p2, = plt.plot(wu_saturation, wu_g_values, 'ro')
	plt.title('g(s) versus saturation')
	plt.xlabel('saturation')
	plt.ylabel('g(s)')
	plt.legend([p1, p2],
			   ["our values", "Wu's values (fitted curve)"], loc='center left', bbox_to_anchor=(1, 0.5), prop = fontP)
	plt.show()
