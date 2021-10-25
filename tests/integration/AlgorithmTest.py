import openpnm as op
import numpy as np


def test_thermal_conduction():
    pass
#    # Generate Network and clean up boundaries (delete z-face pores)
#    divs = [10, 50]
#    Lc = 1/divs[0]
#    pn = op.network.Cubic(shape=divs, spacing=Lc)
#    pn.add_boundary_pores()
#    op.topotools.trim(network=pn,
#                      pores=pn.pores(['top_boundary', 'bottom_boundary']))
#    # Generate Geometry objects for internal and boundary pores
#    Ps = pn.pores('internal')
#    Ts = pn.throats()
#    geom = op.geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)
#    geom['pore.area'] = Lc**2
#    geom['pore.diameter'] = Lc
#    geom['throat.length'] = 1e-25
#    geom['throat.area'] = Lc**2
#    Ps = pn.pores('*boundary')
#    boun = op.geometry.GenericGeometry(network=pn, pores=Ps)
#    boun['pore.area'] = Lc**2
#    boun['pore.diameter'] = 1e-25
#    # Create Phase object and associate with a Physics object
#    Cu = op.phases.GenericPhase(network=pn)
#    Cu['pore.thermal_conductivity'] = 1.0  # W/m.K
#    phys1 = op.physics.GenericPhysics(network=pn, phase=Cu, geometry=geom)
#    phys2 = op.physics.GenericPhysics(network=pn, phase=Cu, geometry=boun)
#    mod = op.models.physics.thermal_conductance.series_resistors
#    phys1.add_model(propname='throat.thermal_conductance', model=mod)
#    phys2.models = phys1.models.copy()
#    phys1.regenerate_models()  # Update the conductance values
#    phys2.regenerate_models()  # Update the conductance values
#    # Setup Algorithm object
#    alg = op.algorithms.FourierConduction(network=pn, phase=Cu)
#    inlets = pn.pores('back_boundary')
#    outlets = pn.pores(['front_boundary', 'left_boundary', 'right_boundary'])
#    T_in = 30*np.sin(np.pi*pn['pore.coords'][inlets, 1]/5)+50
#    alg.set_value_BC(values=T_in, pores=inlets)
#    alg.set_value_BC(values=50, pores=outlets)
#    alg.run()
#    Cu.update(alg.results())
#    # Calculate analytical solution over the same domain spacing
#    T = 30*np.sinh(np.pi*pn['pore.coords'][:, 0]/5)/np.sinh(np.pi/5) * \
#        np.sin(np.pi*pn['pore.coords'][:, 1]/5) + 50
#    Cu['pore.analytical_temp'] = T
#    b = Cu['pore.analytical_temp'][pn.pores(geom.name)]
#    a = Cu['pore.temperature'][pn.pores(geom.name)]
#    a = np.reshape(a, (divs[0], divs[1]))
#    b = np.reshape(b, (divs[0], divs[1]))
#    diff = a - b
#    assert np.amax(np.absolute(diff)) < 0.015


def test_open_air_diffusivity():
    pass
#    pn = op.network.Cubic([5, 5, 5], spacing=1)
#    pn.add_boundary_pores()
#    air = op.phases.Air(network=pn)
#    Dab = np.mean(air['pore.diffusivity'])
#    c = np.mean(air['pore.molar_density'])
#    air['throat.diffusive_conductance'] = Dab * c
#    BC1_pores = pn.pores(labels=['top_boundary'])
#    BC2_pores = pn.pores(labels=['bottom_boundary'])
#    Diff = op.algorithms.FickianDiffusion(network=pn, phase=air)
#    # Assign value boundary conditions to top and bottom surface pores
#    Diff.set_value_BC(values=0.6, pores=BC1_pores)
#    Diff.set_value_BC(values=0.4, pores=BC2_pores)
#    Diff.run()
#    Diff.domain_area = 25
#    Diff.domain_length = 5
#    Diff_deff = Diff.calc_effective_diffusivity()/Dab
#    assert np.round(Diff_deff, 3) == 1


def test_Darcy_alg():
    pass
#    # Generate Network and clean up some of boundaries
#    divs = [1, 50, 10]
#    Lc = 0.00004
#    pn = op.network.Cubic(shape=divs, spacing=Lc)
#    # Generate Geometry objects for internal and boundary pores
#    geom = op.geometry._StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
#    # Create Phase object and associate with a Physics object
#    air = op.phases.Air(network=pn)
#    phys = op.physics.GenericPhysics(network=pn, phase=air, geometry=geom)
#    mod = op.models.physics.hydraulic_conductance.hagen_poiseuille
#    phys.add_model(propname='throat.hydraulic_conductance',
#                   model=mod,
#                   calc_pore_len=False)
#    phys.regenerate_models()  # Update the conductance values
#    # Setup Algorithm objects
#    alg1 = op.algorithms.StokesFlow(network=pn, phase=air)
#    inlets = pn.pores('bottom')
#    outlets = pn.pores('top')
#    P_out = 0  # Pa
#    Q_in = 0.6667*(Lc**2)*divs[1]*divs[0]  # m^3/s
#    alg1.set_rate_BC(values=-Q_in, pores=inlets)
#    alg1.set_value_BC(values=P_out, pores=outlets)
#    alg1.run()
#    air.update(alg1.results())
#    a = round(np.absolute(alg1.rate(outlets))[0], 16)
#    b = round(np.absolute(np.sum(alg1['pore.bc_rate'][inlets])), 16)
#    assert a == b
#
#    alg2 = op.algorithms.StokesFlow(network=pn, phase=air)
#    inlets = pn.pores('bottom')
#    outlets = pn.pores('top')
#    P_out = 0  # Pa
#    P_in = 1000  # Pa
#    alg2.set_value_BC(values=P_in, pores=inlets)
#    alg2.set_value_BC(values=P_out, pores=outlets)
#    alg2.run()
#    a = round(np.absolute(alg2.rate(inlets))[0], 16)
#    b = round(np.absolute(alg2.rate(outlets))[0], 16)
#    assert a == b
#    Q = -alg2.rate(inlets)
#    K = Q*air['pore.viscosity'][0]*divs[2]*Lc/(divs[0]*divs[1]*Lc**2*(P_in-P_out))
#    K_alg = alg2.calc_eff_permeability()
