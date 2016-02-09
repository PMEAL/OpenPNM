import scipy as sp
import OpenPNM
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
        flag = False
        try:
            self.alg.set_inlets(pores=self.net.pores('top'), mode='add')
        except:
            flag = True
        assert flag

    def test_set_outlets_conflicting_with_inlets(self):
        self.alg.setup(invading_phase=self.water,
                       defending_phase=self.air,
                       trapping=True)
        self.alg['pore.inlets'][self.net.pores('top')] = True
        flag = False
        try:
            self.alg.set_outlets(pores=self.net.pores('top'))
        except:
            flag = True
        assert flag

    def test_set_outlets_without_trapping(self):
        self.alg.setup(invading_phase=self.water,
                       defending_phase=self.air)
        flag = False
        try:
            self.alg.set_outlets(pores=self.net.pores('top'))
        except:
            flag = True
        assert flag

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

        flag = False
        try:
            self.alg.set_boundary_conditions(pores=Ps, bc_type='bad_type')
        except:
            flag = True
        assert flag

        flag = False
        try:
            self.alg.set_boundary_conditions(bc_type=None, mode='bad_type')
        except:
            flag = True
        assert flag

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
        flag = False
        try:
            self.alg.run()
        except:
            flag = True
        assert flag

    def test_run_w_trapping_but_no_outlets(self):
        self.alg.setup(invading_phase=self.water,
                       defending_phase=self.air,
                       trapping=True)
        Ps = sp.random.randint(0, self.net.Np, 10)
        self.alg.set_boundary_conditions(pores=Ps, bc_type='inlets')
        flag = False
        try:
            self.alg.run()
        except:
            flag = True
        assert flag

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
