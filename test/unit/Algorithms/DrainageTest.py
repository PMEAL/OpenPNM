import scipy as sp
import OpenPNM
ctrl = OpenPNM.Base.Controller()
ctrl.loglevel = 60


class DrainageTest:
    def setup_test(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        self.geo = OpenPNM.Geometry.Toray090(network=self.net,
                                             pores=self.net.Ps,
                                             throats=self.net.Ts)
        self.phase = OpenPNM.Phases.Water(network=self.net)
        self.phys = OpenPNM.Physics.Standard(network=self.net,
                                             phase=self.phase,
                                             pores=self.net.Ps,
                                             throats=self.net.Ts)
        self.alg = OpenPNM.Algorithms.Drainage(network=self.net)

    def test_set_inlets_modes(self):
        self.alg.setup(invading_phase=self.phase)

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
        self.alg.setup(invading_phase=self.phase)
        self.alg['pore.outlets'][self.net.pores('top')] = True
        flag = False
        try:
            self.alg.set_inlets(pores=self.net.pores('top'), mode='add')
        except:
            flag = True
        assert flag

    def test_set_outlets_modes(self):
        self.alg.setup(invading_phase=self.phase)

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
        self.alg.setup(invading_phase=self.phase)

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
        pass
