import scipy as sp
import matplotlib
import OpenPNM
import pytest
from OpenPNM.Algorithms.__OrdinaryPercolation__ import OrdinaryPercolation
mgr = OpenPNM.Base.Workspace()


def test_IP_old_approach():
    mgr.clear()
    pn = OpenPNM.Network.Cubic(shape=[30, 30, 1], spacing=0.01)
    geom = OpenPNM.Geometry.Toray090(network=pn, pores=pn.Ps, throats=pn.Ts)
    water = OpenPNM.Phases.Water(network=pn)
    phys = OpenPNM.Physics.Standard(network=pn, phase=water, geometry=geom)
    inlets = pn.pores('left')
    IP_1 = OpenPNM.Algorithms.InvasionPercolation(network=pn)
    IP_1.run(phase=water, inlets=inlets)
    a = ['pore.all', 'pore.invaded', 'pore.invasion_sequence', 'throat.all',
         'throat.entry_pressure', 'throat.invaded', 'throat.invasion_sequence',
         'throat.order', 'throat.sorted']
    assert sorted(list(IP_1.keys())) == a
    IP_1.return_results()
    assert 'throat.invasion_sequence' in water.keys()
    assert 'pore.invasion_sequence' in water.keys()
    mgr.clear()


def test_IP_new_approach():
    mgr.clear()
    pn = OpenPNM.Network.Cubic(shape=[30, 30, 1], spacing=0.01)
    geom = OpenPNM.Geometry.Toray090(network=pn, pores=pn.Ps, throats=pn.Ts)
    water = OpenPNM.Phases.Water(network=pn)
    phys = OpenPNM.Physics.Standard(network=pn, phase=water, geometry=geom)
    inlets = pn.pores('left')
    IP_1 = OpenPNM.Algorithms.InvasionPercolation(network=pn)
    IP_1.setup(phase=water)
    IP_1.set_inlets(pores=inlets)
    IP_1.run()
    a = ['pore.all', 'pore.invaded', 'pore.invasion_sequence', 'throat.all',
         'throat.entry_pressure', 'throat.invaded', 'throat.invasion_sequence',
         'throat.order', 'throat.sorted']
    assert sorted(list(IP_1.keys())) == a
    IP_1.return_results()
    assert 'throat.invasion_sequence' in water.keys()
    assert 'pore.invasion_sequence' in water.keys()
    mgr.clear()


def test_OP_old_approach():
    mgr.clear()
    pn = OpenPNM.Network.Cubic(shape=[30, 30, 1], spacing=0.01)
    geom = OpenPNM.Geometry.Toray090(network=pn, pores=pn.Ps, throats=pn.Ts)
    water = OpenPNM.Phases.Water(network=pn)
    phys = OpenPNM.Physics.Standard(network=pn, phase=water, geometry=geom)
    OP_1 = OrdinaryPercolation(network=pn, invading_phase=water)
    Ps = pn.pores(labels=['left'])
    OP_1.run(inlets=Ps)
    OP_1.return_results(Pc=7000)
    a = ['pore.all', 'pore.inlets', 'pore.inv_Pc', 'pore.inv_sat',
         'pore.inv_seq', 'throat.all', 'throat.entry_pressure',
         'throat.inv_Pc', 'throat.inv_sat', 'throat.inv_seq']
    assert sorted(list(OP_1.keys())) == a
    mgr.clear()


def test_OP_new_approach():
    mgr.clear()
    pn = OpenPNM.Network.Cubic(shape=[30, 30, 1], spacing=0.01)
    geom = OpenPNM.Geometry.Toray090(network=pn, pores=pn.Ps, throats=pn.Ts)
    water = OpenPNM.Phases.Water(network=pn)
    phys = OpenPNM.Physics.Standard(network=pn, phase=water, geometry=geom)
    inlets = pn.pores('left')
    OP = OrdinaryPercolation(network=pn)
    OP.setup(invading_phase=water)
    OP.set_inlets(pores=inlets)
    OP.run(npts=25)
    a = ['pore.all', 'pore.inlets', 'pore.inv_Pc', 'pore.inv_sat',
         'pore.inv_seq', 'throat.all', 'throat.entry_pressure',
         'throat.inv_Pc', 'throat.inv_sat', 'throat.inv_seq']
    assert sorted(list(OP.keys())) == a
    V_inv = sp.sum(pn['pore.volume'][OP['pore.inv_Pc'] < sp.inf])
    V_tot = sp.sum(pn['pore.volume'])
    assert V_inv/V_tot == 1.0
    mgr.clear()


def test_OP_trapping():
    mgr.clear()
    pn = OpenPNM.Network.Cubic(shape=[30, 30, 1], spacing=0.01)
    geom = OpenPNM.Geometry.Toray090(network=pn, pores=pn.Ps, throats=pn.Ts)
    water = OpenPNM.Phases.Water(network=pn)
    air = OpenPNM.Phases.Air(network=pn)
    phys = OpenPNM.Physics.GenericPhysics(network=pn,
                                          geometry=geom,
                                          phase=water)
    f = OpenPNM.Physics.models.capillary_pressure.washburn
    phys.models.add(propname='throat.capillary_pressure',
                    model=f)

    OP = OrdinaryPercolation(network=pn,
                             invading_phase=water,
                             defending_phase=air)
    OP.run(inlets=pn.pores('left'), outlets=pn.pores('right'), trapping=True)
    V_inv = sp.sum(pn['pore.volume'][OP['pore.inv_Pc'] < sp.inf])
    V_tot = sp.sum(pn['pore.volume'])
    assert V_inv/V_tot < 1.0
    mgr.clear()


def test_OP_plotting():
    mgr.clear()
    pn = OpenPNM.Network.Cubic(shape=[30, 30, 1], spacing=0.01)
    geom = OpenPNM.Geometry.Toray090(network=pn, pores=pn.Ps, throats=pn.Ts)
    water = OpenPNM.Phases.Water(network=pn)
    OpenPNM.Physics.Standard(network=pn, phase=water, geometry=geom)
    inlets = pn.pores('left')
    OP = OrdinaryPercolation(network=pn)
    OP.setup(invading_phase=water)
    OP.set_inlets(pores=inlets)
    OP.run(npts=25)
    a = OP.plot_drainage_curve()
    assert isinstance(a, matplotlib.figure.Figure)
    a = OP.plot_primary_drainage_curve()
    assert isinstance(a, matplotlib.figure.Figure)


def test_MixedPercolation():
    from OpenPNM.Physics import models as pm
    from OpenPNM.Geometry import models as gm
    import numpy as np
    mgr.clear()
    mgr.loglevel = 50
    fiber_rad = 2e-6
    pn = OpenPNM.Network.Cubic(shape=[5, 5, 5], spacing=5e-5, name='net')
    pn.add_boundaries()
    Ps = pn.pores()
    Ts = pn.throats()
    geom = OpenPNM.Geometry.Toray090(network=pn, pores=Ps, throats=Ts,
                                     name='geom')
    geom.models.add(propname='throat.centroid',
                    model=gm.throat_centroid.pore_coords)
    geom.models.add(propname='throat.normal',
                    model=gm.throat_normal.pore_coords)
    # Phases
    air = OpenPNM.Phases.Air(network=pn, name='air')
    water = OpenPNM.Phases.Water(network=pn, name='water')
    air['pore.contact_angle'] = 70
    water['pore.contact_angle'] = 110
    air["pore.surface_tension"] = water["pore.surface_tension"]
    # Physics
    Ps = pn.pores()
    Ts = pn.throats()
    phys_water = OpenPNM.Physics.Standard(network=pn, phase=water,
                                          pores=Ps, throats=Ts)
    phys_air = OpenPNM.Physics.Standard(network=pn, phase=air,
                                        pores=Ps, throats=Ts)
    throat_diam = 'throat.diameter'
    pore_diam = 'pore.diameter'
    pmod = pm.capillary_pressure.purcell_bi
    phys_air.models.add(propname='throat.capillary_pressure',
                        model=pmod,
                        r_toroid=fiber_rad,
                        diameter=throat_diam,
                        h_max=pore_diam)
    phys_water.models.add(propname='throat.capillary_pressure',
                          model=pmod,
                          r_toroid=fiber_rad,
                          diameter=throat_diam,
                          h_max=pore_diam)
    phys_air.models.add(propname='throat.snap_off',
                        model=pm.capillary_pressure.ransohoff_snap_off,
                        throat_diameter=throat_diam,
                        wavelength=fiber_rad)
    phys_water.models.add(propname='throat.snap_off',
                          model=pm.capillary_pressure.ransohoff_snap_off,
                          throat_diameter=throat_diam,
                          wavelength=fiber_rad)
    phys_air['throat.snap_off'] = np.abs(phys_air['throat.snap_off'])
    phys_air['pore.capillary_pressure'] = 0
    phys_water['pore.capillary_pressure'] = 0
    BPs = pn.pores('boundary')
    NBPs = pn.find_neighbor_pores(BPs, flatten=False)
    boundary_neighbors = []
    for NBP in NBPs:
        boundary_neighbors.append(NBP[0])
    NBPs = np.asarray(boundary_neighbors)
    wPc_NBPs = phys_water["pore.capillary_pressure"][NBPs]
    phys_water["pore.capillary_pressure"][BPs] = wPc_NBPs
    aPc_NBPs = phys_air["pore.capillary_pressure"][NBPs]
    phys_air["pore.capillary_pressure"][BPs] = aPc_NBPs
    # Water Invasion
    inv_points = np.linspace(0, 20000, 11)
    inlets = pn.pores(labels=['bottom_boundary'])
    outlets = pn.pores(labels=['top_boundary'])
    inlet_inv_seq = -1
    IP_1 = OpenPNM.Algorithms.MixedPercolation(network=pn, name='invasion')
    IP_1.setup(phase=water,
               def_phase=air,
               inlets=inlets,
               inlet_inv_seq=inlet_inv_seq,
               snap_off=False,
               coop_fill=True,
               partial=False,
               fill_radius=fiber_rad)
    IP_1.run(inv_points=inv_points)
    IP_1.apply_trapping(outlets=outlets)
    inv_data = IP_1.plot_drainage_curve(inv_points=inv_points, lpf=True)
    # Air Invasion
    IP_1.return_results()
    sat = inv_data[1]
    assert len(sat) == 11
    assert isinstance(inv_data[0], matplotlib.figure.Figure)
