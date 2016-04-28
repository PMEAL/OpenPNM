import scipy as sp
import OpenPNM
import pytest
import matplotlib as mpl
mgr = OpenPNM.Base.Workspace()
mgr.loglevel = 60


class DrainageTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        self.geo = OpenPNM.Geometry.Toray090(network=self.net,
                                             pores=self.net.Ps,
                                             throats=self.net.Ts)
        self.water = OpenPNM.Phases.Water(network=self.net)
        self.air = OpenPNM.Phases.Air(network=self.net)
        self.phys = OpenPNM.Physics.Standard(network=self.net,
                                             phase=self.water,
                                             pores=self.net.Ps,
                                             throats=self.net.Ts)
        self.alg = OpenPNM.Algorithms.Drainage(network=self.net)

    def test_set_inlets_modes(self):
        self.alg.setup(invading_phase=self.water, defending_phase=self.air)

        self.alg.set_inlets(pores=self.net.pores('top'), mode='add')
        assert sp.sum(self.alg['pore.inlets']) == 25

        self.alg.set_inlets(pores=self.net.pores('bottom'), mode='add')
        assert sp.sum(self.alg['pore.inlets']) == 50

        self.alg.set_inlets(pores=self.net.pores('top'), mode='overwrite')
        assert sp.sum(self.alg['pore.inlets']) == 25

        self.alg.set_inlets(pores=self.net.pores('top'), mode='remove')
        assert sp.sum(self.alg['pore.inlets']) == 0

        self.alg.set_inlets(pores=self.net.pores('top'), mode='add')
        self.alg.set_inlets(mode='clear')
        assert sp.sum(self.alg['pore.inlets']) == 0

    def test_set_inlets_conflicting_with_outlets(self):
        self.alg.setup(invading_phase=self.water, defending_phase=self.air)
        self.alg['pore.outlets'][self.net.pores('top')] = True
        with pytest.raises(Exception):
            self.alg.set_inlets(pores=self.net.pores('top'), mode='add')

    def test_set_outlets_conflicting_with_inlets(self):
        self.alg.setup(invading_phase=self.water,
                       defending_phase=self.air,
                       trapping=True)
        self.alg['pore.inlets'][self.net.pores('top')] = True
        with pytest.raises(Exception):
            self.alg.set_outlets(pores=self.net.pores('top'))

    def test_set_outlets_without_trapping(self):
        self.alg.setup(invading_phase=self.water,
                       defending_phase=self.air)
        with pytest.raises(Exception):
            self.alg.set_outlets(pores=self.net.pores('top'))

    def test_set_outlets_modes(self):
        self.alg.setup(invading_phase=self.water,
                       defending_phase=self.air,
                       trapping=True)

        self.alg.set_outlets(pores=self.net.pores('top'), mode='add')
        assert sp.sum(self.alg['pore.outlets']) == 25

        self.alg.set_outlets(pores=self.net.pores('bottom'), mode='add')
        assert sp.sum(self.alg['pore.outlets']) == 50

        self.alg.set_outlets(pores=self.net.pores('top'), mode='overwrite')
        assert sp.sum(self.alg['pore.outlets']) == 25

        self.alg.set_outlets(pores=self.net.pores('top'), mode='remove')
        assert sp.sum(self.alg['pore.outlets']) == 0

        self.alg.set_outlets(pores=self.net.pores('top'), mode='add')
        self.alg.set_outlets(mode='clear')
        assert sp.sum(self.alg['pore.outlets']) == 0

    def test_set_residual_modes(self):
        self.alg.setup(invading_phase=self.water, defending_phase=self.air)

        Ps = sp.random.randint(0, self.net.Np, 10)
        Ts = self.net.find_neighbor_pores(pores=Ps)
        self.alg.set_residual(pores=Ps, throats=Ts, mode='add')
        assert sp.sum(self.alg['pore.residual']) == sp.size(sp.unique(Ps))
        assert sp.sum(self.alg['throat.residual']) == sp.size(sp.unique(Ts))

        Ps = sp.random.randint(0, self.net.Np, 10)
        Ts = self.net.find_neighbor_pores(pores=Ps)
        self.alg.set_residual(pores=Ps, throats=Ts, mode='add')
        assert sp.sum(self.alg['pore.residual']) > sp.size(sp.unique(Ps))
        assert sp.sum(self.alg['throat.residual']) > sp.size(sp.unique(Ts))

        Ps = sp.random.randint(0, self.net.Np, 10)
        Ts = self.net.find_neighbor_pores(pores=Ps)
        self.alg.set_residual(pores=Ps, throats=Ts, mode='overwrite')
        assert sp.sum(self.alg['pore.residual']) == sp.size(sp.unique(Ps))
        assert sp.sum(self.alg['throat.residual']) == sp.size(sp.unique(Ts))

        self.alg.set_residual(pores=Ps, throats=Ts, mode='remove')
        assert sp.sum(self.alg['pore.residual']) == 0

        self.alg.set_residual(pores=Ps, throats=Ts, mode='add')
        self.alg.set_residual(mode='clear')
        assert sp.sum(self.alg['pore.residual']) == 0

    def test_set_boundary_conditions_clear(self):
        self.alg.setup(invading_phase=self.water,
                       defending_phase=self.air,
                       trapping=True)

        self.alg['pore.inlets'] = True
        self.alg['pore.outlets'] = True
        self.alg['pore.residual'] = True
        self.alg['throat.residual'] = True

        self.alg.set_boundary_conditions(mode='clear')
        assert sp.sum(self.alg['pore.inlets']) == 0
        assert sp.sum(self.alg['pore.outlets']) == 0
        assert sp.sum(self.alg['pore.residual']) == 0
        assert sp.sum(self.alg['throat.residual']) == 0

    def test_set_boundary_conditions_bctypes(self):
        self.alg.setup(invading_phase=self.water,
                       defending_phase=self.air,
                       trapping=True)
        Ps = sp.random.randint(0, self.net.Np, 10)

        self.alg.set_boundary_conditions(pores=Ps, bc_type='inlets')
        assert sp.sum(self.alg['pore.inlets']) == sp.size(sp.unique(Ps))
        self.alg['pore.inlets'] = False

        self.alg.set_boundary_conditions(pores=Ps, bc_type='outlets')
        assert sp.sum(self.alg['pore.outlets']) == sp.size(sp.unique(Ps))
        self.alg['pore.outlets'] = False

        self.alg.set_boundary_conditions(pores=Ps, bc_type='residual')
        assert sp.sum(self.alg['pore.residual']) == sp.size(sp.unique(Ps))
        self.alg['pore.residual'] = False

        with pytest.raises(Exception):
            self.alg.set_boundary_conditions(pores=Ps, bc_type='bad_type')

        with pytest.raises(Exception):
            self.alg.set_boundary_conditions(bc_type=None, mode='bad_type')

    def test_run_npts(self):
        self.alg.setup(invading_phase=self.water, defending_phase=self.air)
        Ps = sp.random.randint(0, self.net.Np, 10)
        self.alg.set_boundary_conditions(pores=Ps, bc_type='inlets')
        self.alg.run(npts=20)

    def test_run_inv_pressures(self):
        self.alg.setup(invading_phase=self.water, defending_phase=self.air)
        Ps = sp.random.randint(0, self.net.Np, 10)
        self.alg.set_boundary_conditions(pores=Ps, bc_type='inlets')
        self.alg.run(inv_pressures=range(0, 20000, 1000))
        assert sp.all(self.alg._inv_points == range(0, 20000, 1000))

    def test_run_no_inlets(self):
        self.alg.setup(invading_phase=self.water, defending_phase=self.air)
        with pytest.raises(Exception):
            self.alg.run()

    def test_run_w_trapping_but_no_outlets(self):
        self.alg.setup(invading_phase=self.water,
                       defending_phase=self.air,
                       trapping=True)
        Ps = sp.random.randint(0, self.net.Np, 10)
        self.alg.set_boundary_conditions(pores=Ps, bc_type='inlets')
        with pytest.raises(Exception):
            self.alg.run()

    def test_run_w_trapping(self):
        self.alg.setup(invading_phase=self.water,
                       defending_phase=self.air,
                       trapping=True)
        self.alg.set_boundary_conditions(pores=self.net.pores('top'),
                                         bc_type='inlets')
        self.alg.set_boundary_conditions(pores=self.net.pores('bottom'),
                                         bc_type='outlets')
        self.alg.run()
        data = self.alg.get_drainage_data()
        assert 'capillary_pressure' in data.keys()
        assert 'invading_phase_saturation' in data.keys()

    def test_run_w_residual_pores_and_throats(self):
        self.alg.setup(invading_phase=self.water, defending_phase=self.air)
        self.alg.set_boundary_conditions(pores=self.net.pores('top'),
                                         bc_type='inlets')
        self.alg.set_boundary_conditions(pores=self.net.pores('bottom'),
                                         bc_type='residual')
        self.alg.run()
        data = self.alg.get_drainage_data()
        assert 'capillary_pressure' in data.keys()
        assert 'invading_phase_saturation' in data.keys()

    def test_basic(self):
        self.alg.setup(invading_phase=self.water, defending_phase=self.air)
        self.alg.set_inlets(pores=self.net.pores('boundary_top'))
        self.alg.run()
        data = self.alg.rainage.get_drainage_data()

        assert sp.amin(data['invading_phase_saturation']) == 0.0
        assert sp.amax(data['invading_phase_saturation']) == 1.0

    def test_residual(self):
        self.alg.setup(invading_phase=self.water, defending_phase=self.air)
        self.alg.set_inlets(pores=self.net.pores('boundary_top'))
        Ps = sp.random.randint(0, self.Np, 1000)
        Ts = sp.random.randint(0, self.Nt, 1000)
        self.alg.set_residual(pores=Ps, throats=Ts)
        self.alg.run()
        data = self.alg.get_drainage_data()

        assert sp.amin(data['invading_phase_saturation']) > 0
        assert sp.amax(data['invading_phase_saturation']) == 1.0

    def test_trapping(self):
        self.alg.setup(invading_phase=self.water, defending_phase=self.air,
                       trapping=True)
        self.alg.set_inlets(pores=self.net.pores('boundary_top'))
        self.alg.set_outlets(pores=self.pores('boundary_bottom')[0:300])
        self.alg.run()
        data = self.alg.get_drainage_data()

        assert sp.amin(data['invading_phase_saturation']) == 0.0
        assert sp.amax(data['invading_phase_saturation']) < 1.0

    def test_late_pore_filling(self):
        mod = OpenPNM.Physics.models.multiphase.late_pore_filling
        self.phys.models.add(propname='pore.fractional_filling',
                             model=mod, Pc=0, Swp_star=0.2, eta=1)
        self.phys.regenerate()
        self.alg.setup(invading_phase=self.water, defending_phase=self.air,
                       pore_filling='pore.fractional_filling')
        self.alg.set_inlets(pores=self.net.pores('boundary_top'))
        self.alg.run()
        data = self.alg.get_drainage_data()
        assert sp.amin(data['invading_phase_saturation']) == 0.0
        assert sp.amax(data['invading_phase_saturation']) < 1.0

        self.alg.return_results(Pc=5000)
        assert 'pore.occupancy' in self.water.keys()
        assert 'pore.partial_occupancy' in self.water.keys()

    def test_late_throat_filling(self):
        mod = OpenPNM.Physics.models.multiphase.late_throat_filling
        self.phys.models.add(propname='throat.fractional_filling',
                             model=mod, Pc=0, Swp_star=0.2, eta=1)
        self.phys.regenerate()
        self.alg.setup(invading_phase=self.water, defending_phase=self.air,
                       throat_filling='throat.fractional_filling')
        self.alg.set_inlets(pores=self.net.pores('boundary_top'))
        self.alg.run()
        data = self.alg.get_drainage_data()

        assert sp.amin(data['invading_phase_saturation']) == 0.0
        assert sp.amax(data['invading_phase_saturation']) < 1.0

        self.alg.return_results(Pc=5000)
        assert 'throat.occupancy' in self.water.keys()
        assert 'throat.partial_occupancy' in self.water.keys()

    def test_late_pore_and_throat_filling(self):
        mod = OpenPNM.Physics.models.multiphase.late_pore_filling
        self.phys.models.add(propname='pore.fractional_filling',
                             model=mod, Pc=0, Swp_star=0.2, eta=1)
        mod = OpenPNM.Physics.models.multiphase.late_throat_filling
        self.phys.models.add(propname='throat.fractional_filling',
                             model=mod, Pc=0, Swp_star=0.2, eta=1)
        self.phys.regenerate()
        self.alg.setup(invading_phase=self.water, defending_phase=self.air,
                       pore_filling='pore.fractional_filling',
                       throat_filling='throat.fractional_filling')
        self.alg.set_inlets(pores=self.net.pores('boundary_top'))
        self.alg.run()
        data = self.alg.get_drainage_data()
        assert sp.amin(data['invading_phase_saturation']) == 0.0
        assert sp.amax(data['invading_phase_saturation']) < 1.0

        self.alg.return_results(Pc=5000)
        assert 'pore.occupancy' in self.water.keys()
        assert 'throat.occupancy' in self.water.keys()
        assert 'pore.partial_occupancy' in self.water.keys()
        assert 'throat.partial_occupancy' in self.water.keys()

    def test_ploting(self):
        self.alg.setup(invading_phase=self.water, defending_phase=self.air)
        self.alg.set_inlets(pores=self.net.pores('boundary_top'))
        self.alg.run()
        data = self.alg.get_drainage_data()
        a = self.alg.plot_drainage_curve(data)
        assert isinstance(a, mpl.figure.Figure)

    def test_residual_and_lpf(self):
        mod = OpenPNM.Physics.models.multiphase.late_pore_filling
        self.phys.models.add(propname='pore.fractional_filling',
                             model=mod, Pc=0, Swp_star=0.2, eta=1)
        mod = OpenPNM.Physics.models.multiphase.late_throat_filling
        self.phys.models.add(propname='throat.fractional_filling',
                             model=mod, Pc=0, Swp_star=0.2, eta=1)
        self.phys.regenerate()
        self.alg.setup(invading_phase=self.water, defending_phase=self.air,
                       pore_filling='pore.fractional_filling',
                       throat_filling='throat.fractional_filling')
        self.alg.set_inlets(pores=self.net.pores('boundary_top'))
        mask = sp.random.random(len(self.net.pores('internal'))) < 0.1
        resPs = self.net.pores('internal')[mask]
        mask = sp.random.random(len(self.net.throats('internal'))) < 0.1
        resTs = self.net.throats('internal')[mask]
        self.alg.set_residual(pores=resPs, throats=resTs)
        self.alg.run()
        self.alg.return_results(Pc=5000)
        data = self.alg.get_drainage_data()
        assert sp.all(self.water["pore.partial_occupancy"][resPs] == 1.0)
        assert sp.all(self.water["throat.partial_occupancy"][resTs] == 1.0)
        assert sp.amin(data['invading_phase_saturation']) > 0.0
        assert sp.amax(data['invading_phase_saturation']) < 1.0
        assert sp.all(self.water["pore.occupancy"] +
                      self.air["pore.occupancy"] == 1.0)
        total_pp = self.water["pore.partial_occupancy"] + self.air["pore.partial_occupancy"]
        assert sp.all(total_pp == 1.0)
        assert sp.all(self.water["throat.occupancy"] +
                      self.air["throat.occupancy"] == 1.0)
        total_pt = self.water["throat.partial_occupancy"] + self.air["throat.partial_occupancy"]
        assert sp.all(total_pt == 1.0)
