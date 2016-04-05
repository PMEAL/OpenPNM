import numpy as np
import scipy as sp
import OpenPNM
import pytest


def test_linear_solvers():
    pn = OpenPNM.Network.Cubic([1, 40, 30], spacing=0.0001)
    geom = OpenPNM.Geometry.Toray090(network=pn,
                                     pores=pn.pores(),
                                     throats=pn.throats())
    air = OpenPNM.Phases.Air(network=pn)
    phys_air = OpenPNM.Physics.Standard(network=pn,
                                        phase=air,
                                        pores=pn.pores(),
                                        throats=pn.throats())

    BC1_pores = pn.pores(labels=['left'])
    BC2_pores = pn.pores(labels=['right'])

    alg_1 = OpenPNM.Algorithms.FickianDiffusion(network=pn, phase=air)
    alg_1.set_boundary_conditions(bctype='Dirichlet',
                                  bcvalue=1,
                                  pores=BC1_pores)
    alg_1.set_boundary_conditions(bctype='Dirichlet',
                                  bcvalue=0,
                                  pores=BC2_pores)
    alg_1.run(iterative_solver='gmres')

    alg_2 = OpenPNM.Algorithms.FickianDiffusion(network=pn, phase=air)
    alg_2.set_boundary_conditions(bctype='Neumann',
                                  bcvalue=-1e-11,
                                  pores=BC1_pores)
    alg_2.set_boundary_conditions(bctype='Dirichlet',
                                  bcvalue=0,
                                  pores=BC2_pores)
    alg_2.run(iterative_solver='cg')

    alg_3 = OpenPNM.Algorithms.FickianDiffusion(network=pn, phase=air)
    alg_3.set_boundary_conditions(bctype='Neumann_group',
                                  bcvalue=-3e-10,
                                  pores=BC1_pores)
    alg_3.set_boundary_conditions(bctype='Dirichlet',
                                  bcvalue=0,
                                  pores=BC2_pores)
    alg_3.run()

    alg_4 = OpenPNM.Algorithms.FickianDiffusion(network=pn, phase=air)
    alg_4.set_boundary_conditions(bctype='Neumann_group',
                                  bcvalue=-3e-10,
                                  pores=BC1_pores)
    alg_4.set_boundary_conditions(bctype='Dirichlet',
                                  bcvalue=0,
                                  pores=BC2_pores)
    alg_4.setup()
    alg_4.solve()

    assert round(sp.absolute(alg_1.rate(BC1_pores))[0], 16) ==\
        round(sp.absolute(alg_1.rate(BC2_pores))[0], 16)
    assert round(sp.absolute(alg_2.rate(BC2_pores))[0], 16) ==\
        round(sp.absolute(sp.unique(alg_2['pore.bcval_Neumann']))[0] *
              len(BC1_pores), 16)
    assert round(sp.absolute(alg_3.rate(BC2_pores))[0], 16) ==\
        round(sp.absolute(sp.unique(alg_3['pore.bcval_Neumann_group']))[0], 16)
    assert round(sp.absolute(alg_4.rate(BC2_pores))[0], 16) ==\
        round(sp.absolute(sp.unique(alg_4['pore.bcval_Neumann_group']))[0], 16)

    assert round(sp.absolute(sp.sum(alg_1.rate(BC1_pores, mode='single'))), 16) ==\
        round(sp.absolute(alg_1.rate(BC1_pores))[0], 16)
    assert round(sp.absolute(sp.sum(alg_2.rate(BC2_pores, mode='single'))), 16) ==\
        round(sp.absolute(alg_2.rate(BC2_pores))[0], 16)
    assert round(sp.absolute(sp.sum(alg_3.rate(BC2_pores, mode='single'))), 16) ==\
        round(sp.absolute(alg_3.rate(BC2_pores))[0], 16)
    assert round(sp.absolute(sp.sum(alg_4.rate(BC2_pores, mode='single'))), 16) ==\
        round(sp.absolute(alg_4.rate(BC2_pores))[0], 16)


def test_add_boundary():
    pn = OpenPNM.Network.Cubic(shape=[5, 5, 5])
    pn.add_boundaries()

    keys_expected = {'pore.back', 'pore.bottom', 'pore.top_boundary',
                     'pore.right_boundary', 'throat.back_boundary',
                     'throat.all', 'throat.bottom_boundary',
                     'throat.front_boundary', 'pore.boundary',
                     'throat.left_boundary', 'throat.conns',
                     'throat.top_boundary', 'pore.back_boundary', 'pore.top',
                     'pore.front_boundary', 'pore.all', 'pore.front',
                     'pore.left_boundary', 'throat.boundary',
                     'pore.bottom_boundary', 'throat.right_boundary',
                     'pore.coords', 'pore.internal', 'pore.index', 'pore.left',
                     'pore.right'}

    keys_found = set(pn.keys())

    symmetric_diff = keys_found ^ keys_expected
    assert not symmetric_diff


def test_open_air_diffusivity():
    pn = OpenPNM.Network.Cubic([5, 5, 5], spacing=1)
    pn.add_boundaries()
    Ps = pn.pores('boundary', mode='not')
    Ts = pn.find_neighbor_throats(pores=Ps, mode='intersection', flatten=True)
    geom = OpenPNM.Geometry.Cube_and_Cuboid(network=pn, pores=Ps, throats=Ts)
    geom['pore.diameter'] = 0.999999
    geom['throat.diameter'] = 0.999999
    geom.regenerate(['pore.diameter', 'throat.diameter'], mode='exclude')
    Ps = pn.pores('boundary')
    Ts = pn.find_neighbor_throats(pores=Ps, mode='not_intersection')
    boun = OpenPNM.Geometry.Boundary(network=pn,
                                     pores=Ps,
                                     throats=Ts,
                                     shape='cubes')
    air = OpenPNM.Phases.Air(network=pn)
    Ps = pn.pores()
    Ts = pn.throats()
    phys_air = OpenPNM.Physics.Standard(network=pn,
                                        phase=air,
                                        pores=Ps,
                                        throats=Ts)
    BC1_pores = pn.pores(labels=['top_boundary'])
    BC2_pores = pn.pores(labels=['bottom_boundary'])
    Diff = OpenPNM.Algorithms.FickianDiffusion(network=pn,
                                               phase=air)
    # Assign Dirichlet boundary conditions to top and bottom surface pores
    Diff.set_boundary_conditions(bctype='Dirichlet',
                                 bcvalue=0.6,
                                 pores=BC1_pores)
    Diff.set_boundary_conditions(bctype='Dirichlet',
                                 bcvalue=0.4,
                                 pores=BC2_pores)
    Diff.run()
    Diff.return_results()
    Diff_deff = Diff.calc_eff_diffusivity()/np.mean(air['pore.diffusivity'])
    assert np.round(Diff_deff, 3) == 1


def test_thermal_conduction():
    # Generate Network and clean up boundaries (delete z-face pores)
    divs = [10, 50]
    Lc = 0.1  # cm
    pn = OpenPNM.Network.Cubic(shape=divs, spacing=Lc)
    pn.add_boundaries()
    pn.trim(pores=pn.pores(['top_boundary', 'bottom_boundary']))
    # Generate Geometry objects for internal and boundary pores
    Ps = pn.pores('internal')
    Ts = pn.throats()
    geom = OpenPNM.Geometry.GenericGeometry(network=pn,
                                            pores=Ps,
                                            throats=Ts)
    geom['pore.area'] = Lc**2
    geom['pore.diameter'] = Lc
    geom['throat.length'] = 1e-25
    geom['throat.area'] = Lc**2
    Ps = pn.pores('boundary')
    boun = OpenPNM.Geometry.GenericGeometry(network=pn, pores=Ps)
    boun['pore.area'] = Lc**2
    boun['pore.diameter'] = 1e-25
    # Create Phase object and associate with a Physics object
    Cu = OpenPNM.Phases.GenericPhase(network=pn)
    Cu['pore.thermal_conductivity'] = 1.0  # W/m.K
    phys = OpenPNM.Physics.GenericPhysics(network=pn,
                                          phase=Cu,
                                          pores=pn.pores(),
                                          throats=pn.throats())
    mod = OpenPNM.Physics.models.thermal_conductance.series_resistors
    phys.add_model(propname='throat.thermal_conductance', model=mod)
    phys.regenerate()  # Update the conductance values
    # Setup Algorithm object
    Fourier_alg = OpenPNM.Algorithms.FourierConduction(network=pn, phase=Cu)
    inlets = pn.pores('back_boundary')
    outlets = pn.pores(['front_boundary', 'left_boundary', 'right_boundary'])
    T_in = 30*sp.sin(sp.pi*pn['pore.coords'][inlets, 1]/5)+50
    Fourier_alg.set_boundary_conditions(bctype='Dirichlet',
                                        bcvalue=T_in,
                                        pores=inlets)
    Fourier_alg.set_boundary_conditions(bctype='Dirichlet',
                                        bcvalue=50,
                                        pores=outlets)
    Fourier_alg.run()
    Fourier_alg.return_results()
    # Calculate analytical solution over the same domain spacing
    Cu['pore.analytical_temp'] = 30*sp.sinh(sp.pi*pn['pore.coords'][:, 0]/5)/sp.sinh(sp.pi/5)*sp.sin(sp.pi*pn['pore.coords'][:, 1]/5) + 50
    b = Cu['pore.analytical_temp'][pn.pores(geom.name)]
    a = Cu['pore.temperature'][pn.pores(geom.name)]
    a = sp.reshape(a, (divs[0], divs[1]))
    b = sp.reshape(b, (divs[0], divs[1]))
    diff = a - b
    assert sp.amax(np.absolute(diff)) < 0.015


def test_Darcy_alg():
    # Generate Network and clean up some of boundaries
    divs = [1, 50, 10]
    Lc = 0.00004
    pn = OpenPNM.Network.Cubic(shape=divs, spacing=Lc)
    pn.add_boundaries()
    Ps = pn.pores(['front_boundary', 'back_boundary'])
    pn.trim(pores=Ps)
    # Generate Geometry objects for internal and boundary pores
    Ps = pn.pores('boundary', mode='not')
    Ts = pn.find_neighbor_throats(pores=Ps, mode='intersection', flatten=True)
    geom = OpenPNM.Geometry.Toray090(network=pn, pores=Ps, throats=Ts)
    Ps = pn.pores('boundary')
    Ts = pn.find_neighbor_throats(pores=Ps, mode='not_intersection')
    boun = OpenPNM.Geometry.Boundary(network=pn, pores=Ps, throats=Ts)
    # Create Phase object and associate with a Physics object
    air = OpenPNM.Phases.Air(network=pn)
    Ps = pn.pores()
    Ts = pn.throats()
    phys = OpenPNM.Physics.GenericPhysics(network=pn,
                                          phase=air,
                                          pores=Ps,
                                          throats=Ts)
    from OpenPNM.Physics import models as pm
    phys.add_model(propname='throat.hydraulic_conductance',
                   model=pm.hydraulic_conductance.hagen_poiseuille,
                   calc_pore_len=False)
    phys.regenerate()  # Update the conductance values
    # Setup Algorithm objects
    Darcy1 = OpenPNM.Algorithms.StokesFlow(network=pn, phase=air)
    inlets = pn.pores('bottom_boundary')
    Ps = pn.pores('top_boundary')
    outlets = Ps[pn['pore.coords'][Ps, 1] < (divs[1]*Lc/2)]
    P_out = 0  # Pa
    Q_in = 0.6667*(Lc**2)*divs[1]*divs[0]  # m^3/s
    Darcy1.set_boundary_conditions(bctype='Neumann_group',
                                   bcvalue=-Q_in,
                                   pores=inlets)
    Darcy1.set_boundary_conditions(bctype='Dirichlet',
                                   bcvalue=P_out,
                                   pores=outlets)
    Darcy1.run()
    Darcy1.return_results()
    print('pore pressure for Darcy1 algorithm:')
    print(air['pore.pressure'])
    Darcy2 = OpenPNM.Algorithms.StokesFlow(network=pn, phase=air)
    inlets = pn.pores('bottom_boundary')
    outlets = pn.pores('top_boundary')
    P_out = 10  # Pa
    P_in = 1000  # Pa
    Darcy2.set_boundary_conditions(bctype='Dirichlet',
                                   bcvalue=P_in,
                                   pores=inlets)
    Darcy2.set_boundary_conditions(bctype='Dirichlet',
                                   bcvalue=P_out,
                                   pores=outlets)
    Darcy2.run()
    print('pore pressure for Darcy2 algorithm:')
    print(Darcy2['pore.pressure'])
    Q = -Darcy2.rate(inlets)
    K = Q*air['pore.viscosity'][0]*divs[2]*Lc/(divs[0]*divs[1]*Lc**2*(P_in-P_out))
    Vp = sp.sum(pn['pore.volume']) + sp.sum(pn['throat.volume'])
    Vb = sp.prod(divs)*Lc**3
    e = Vp/Vb
    print('Effective permeability: ', K, '- Porosity: ', e)

    a = round(sp.absolute(Darcy1.rate(outlets))[0], 16)
    pore_prop = 'pore.bcval_Neumann_group'
    b = round(sp.absolute(sp.unique(Darcy1[pore_prop]))[0], 16)
    assert a == b

    a = round(sp.absolute(Darcy2.rate(inlets))[0], 16)
    b = round(sp.absolute(Darcy2.rate(outlets))[0], 16)
    assert a == b
