import scipy as sp
import OpenPNM
import matplotlib

pn = OpenPNM.Network.Cubic(shape=[25, 25, 25], spacing=0.00005)
pn['pore.internal'] = True
pn['throat.internal'] = True
pn.add_boundary_pores(pores=pn.pores('top'),
                      offset=[0, 0, 0.0001],
                      apply_label='boundary_top')
pn.add_boundary_pores(pores=pn.pores('bottom'),
                      offset=[0, 0, -0.0001],
                      apply_label='boundary_bottom')
geom = OpenPNM.Geometry.Toray090(network=pn,
                                 pores=pn.pores('internal'),
                                 throats=pn.throats('internal'))
boun = OpenPNM.Geometry.Boundary(network=pn,
                                 pores=pn.pores('internal', mode='not'),
                                 throats=pn.throats('internal', mode='not'))
water = OpenPNM.Phases.Water(network=pn)
air = OpenPNM.Phases.Air(network=pn)
phys = OpenPNM.Physics.Standard(network=pn,
                                phase=water,
                                pores=pn.Ps,
                                throats=pn.Ts)
drainage = OpenPNM.Algorithms.Drainage(network=pn)


# def test_basic():
drainage.setup(invading_phase=water, defending_phase=air)
drainage.set_inlets(pores=pn.pores('boundary_top'))
drainage.run()
data = drainage.get_drainage_data()

assert sp.amin(data['invading_phase_saturation']) == 0.0
assert sp.amax(data['invading_phase_saturation']) == 1.0


# def test_residual():
drainage.setup(invading_phase=water, defending_phase=air)
drainage.set_inlets(pores=pn.pores('boundary_top'))
Ps = sp.random.randint(0, pn.Np, 1000)
Ts = sp.random.randint(0, pn.Nt, 1000)
drainage.set_residual(pores=Ps, throats=Ts)
drainage.run()
data = drainage.get_drainage_data()

assert sp.amin(data['invading_phase_saturation']) > 0
assert sp.amax(data['invading_phase_saturation']) == 1.0


# def test_trapping():
drainage.setup(invading_phase=water, defending_phase=air, trapping=True)
drainage.set_inlets(pores=pn.pores('boundary_top'))
drainage.set_outlets(pores=pn.pores('boundary_bottom')[0:300])
drainage.run()
data = drainage.get_drainage_data()

assert sp.amin(data['invading_phase_saturation']) == 0.0
assert sp.amax(data['invading_phase_saturation']) < 1.0


# def test_late_pore_filling():
phys.models.add(propname='pore.fractional_filling',
                model=OpenPNM.Physics.models.multiphase.late_pore_filling,
                Pc=0,
                Swp_star=0.2,
                eta=1)
phys.regenerate()
drainage.setup(invading_phase=water, defending_phase=air,
               pore_filling='pore.fractional_filling')
drainage.set_inlets(pores=pn.pores('boundary_top'))
drainage.run()
data = drainage.get_drainage_data()
assert sp.amin(data['invading_phase_saturation']) == 0.0
assert sp.amax(data['invading_phase_saturation']) < 1.0

drainage.return_results(Pc=5000)
assert 'pore.occupancy' in water.keys()
assert 'pore.partial_occupancy' in water.keys()


# def test_late_throat_filling():
mod = OpenPNM.Physics.models.multiphase.late_throat_filling
phys.models.add(propname='throat.fractional_filling',
                model=mod,
                Pc=0,
                Swp_star=0.2,
                eta=1)
phys.regenerate()
drainage.setup(invading_phase=water, defending_phase=air,
               throat_filling='throat.fractional_filling')
drainage.set_inlets(pores=pn.pores('boundary_top'))
drainage.run()
data = drainage.get_drainage_data()

assert sp.amin(data['invading_phase_saturation']) == 0.0
assert sp.amax(data['invading_phase_saturation']) < 1.0

drainage.return_results(Pc=5000)
assert 'throat.occupancy' in water.keys()
assert 'throat.partial_occupancy' in water.keys()


# def test_late_pore_and_throat_filling():
phys.models.add(propname='pore.fractional_filling',
                model=OpenPNM.Physics.models.multiphase.late_pore_filling,
                Pc=0,
                Swp_star=0.2,
                eta=1)
mod = OpenPNM.Physics.models.multiphase.late_throat_filling
phys.models.add(propname='throat.fractional_filling',
                model=mod,
                Pc=0,
                Swp_star=0.2,
                eta=1)
phys.regenerate()
drainage.setup(invading_phase=water, defending_phase=air,
               pore_filling='pore.fractional_filling',
               throat_filling='throat.fractional_filling')
drainage.set_inlets(pores=pn.pores('boundary_top'))
drainage.run()
data = drainage.get_drainage_data()
assert sp.amin(data['invading_phase_saturation']) == 0.0
assert sp.amax(data['invading_phase_saturation']) < 1.0

drainage.return_results(Pc=5000)
assert 'pore.occupancy' in water.keys()
assert 'throat.occupancy' in water.keys()
assert 'pore.partial_occupancy' in water.keys()
assert 'throat.partial_occupancy' in water.keys()


# def test_ploting():
drainage.setup(invading_phase=water, defending_phase=air)
drainage.set_inlets(pores=pn.pores('boundary_top'))
drainage.run()
data = drainage.get_drainage_data()
a = drainage.plot_drainage_curve(data)
assert isinstance(a, matplotlib.figure.Figure)


# def test_residual_and_lpf():
phys.models.add(propname='pore.fractional_filling',
                model=OpenPNM.Physics.models.multiphase.late_pore_filling,
                Pc=0,
                Swp_star=0.2,
                eta=1)
phys.models.add(propname='throat.fractional_filling',
                model=OpenPNM.Physics.models.multiphase.late_throat_filling,
                Pc=0,
                Swp_star=0.2,
                eta=1)
phys.regenerate()
drainage.setup(invading_phase=water, defending_phase=air,
               pore_filling='pore.fractional_filling',
               throat_filling='throat.fractional_filling')
drainage.set_inlets(pores=pn.pores('boundary_top'))
mask = sp.random.random(len(pn.pores('internal'))) < 0.1
resPs = pn.pores('internal')[mask]
mask = sp.random.random(len(pn.throats('internal'))) < 0.1
resTs = pn.throats('internal')[mask]
drainage.set_residual(pores=resPs, throats=resTs)
drainage.run()
drainage.return_results(Pc=5000)
data = drainage.get_drainage_data()
assert sp.all(water["pore.partial_occupancy"][resPs] == 1.0)
assert sp.all(water["throat.partial_occupancy"][resTs] == 1.0)
assert sp.amin(data['invading_phase_saturation']) > 0.0
assert sp.amax(data['invading_phase_saturation']) < 1.0
assert sp.all(water["pore.occupancy"]+air["pore.occupancy"] == 1.0)
total_pp = water["pore.partial_occupancy"]+air["pore.partial_occupancy"]
assert sp.all(total_pp == 1.0)
assert sp.all(water["throat.occupancy"]+air["throat.occupancy"] == 1.0)
total_pt = water["throat.partial_occupancy"]+air["throat.partial_occupancy"]
assert sp.all(total_pt == 1.0)
