import openpnm as op
import numpy as np


def test_thermal_conduction():
    # Generate Network and clean up boundaries (delete z-face pores)
    divs = [10, 50]
    Lc = 0.1  # cm
    pn = op.network.Cubic(shape=divs, spacing=Lc)
    pn.add_boundary_pores()
    op.topotools.trim(network=pn,
                      pores=pn.pores(['top_boundary', 'bottom_boundary']))
    # Generate Geometry objects for internal and boundary pores
    Ps = pn.pores('internal')
    Ts = pn.throats()
    geom = op.geometry.GenericGeometry(network=pn, pores=Ps, throats=Ts)
    geom['pore.area'] = Lc**2
    geom['pore.diameter'] = Lc
    geom['throat.length'] = 1e-25
    geom['throat.area'] = Lc**2
    Ps = pn.pores('*boundary')
    boun = op.geometry.GenericGeometry(network=pn, pores=Ps)
    boun['pore.area'] = Lc**2
    boun['pore.diameter'] = 1e-25
    # Create Phase object and associate with a Physics object
    Cu = op.phases.GenericPhase(network=pn)
    Cu['pore.thermal_conductivity'] = 1.0  # W/m.K
    phys1 = op.physics.GenericPhysics(network=pn, phase=Cu, geometry=geom)
    phys2 = op.physics.GenericPhysics(network=pn, phase=Cu, geometry=boun)
    mod = op.models.physics.thermal_conductance.series_resistors
    phys1.add_model(propname='throat.thermal_conductance', model=mod)
    phys2.models = phys1.models.copy()
    phys1.regenerate_models()  # Update the conductance values
    phys2.regenerate_models()  # Update the conductance values
    # Setup Algorithm object
    alg = op.algorithms.FourierConduction(network=pn, phase=Cu)
    inlets = pn.pores('back_boundary')
    outlets = pn.pores(['front_boundary', 'left_boundary', 'right_boundary'])
    T_in = 30*np.sin(np.pi*pn['pore.coords'][inlets, 1]/5)+50
    alg.set_dirichlet_BC(values=T_in, pores=inlets)
    alg.set_dirichlet_BC(values=50, pores=outlets)
    alg.run()
    Cu.update(alg.results())
    # Calculate analytical solution over the same domain spacing
    T = 30*np.sinh(np.pi*pn['pore.coords'][:, 0]/5)/np.sinh(np.pi/5) * \
        np.sin(np.pi*pn['pore.coords'][:, 1]/5) + 50
    Cu['pore.analytical_temp'] = T
    b = Cu['pore.analytical_temp'][pn.pores(geom.name)]
    a = Cu['pore.temperature'][pn.pores(geom.name)]
    a = np.reshape(a, (divs[0], divs[1]))
    b = np.reshape(b, (divs[0], divs[1]))
    diff = a - b
    assert np.amax(np.absolute(diff)) < 0.015
