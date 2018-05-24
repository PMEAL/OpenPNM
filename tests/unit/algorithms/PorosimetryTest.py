import openpnm as op
import scipy as sp
import pytest
mgr = op.Workspace()


class PorosimetryTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=0.0005)
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)
        self.water = op.phases.Water(network=self.net)
        self.air = op.phases.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.water,
                                              geometry=self.geo)
        mod = op.models.physics.capillary_pressure.washburn
        self.phys.add_model(propname='throat.entry_pressure',
                            model=mod)

    def test_set_inlets_overwrite(self):
        self.alg = op.algorithms.Porosimetry(network=self.net)
        self.alg.setup(phase=self.water)

        self.alg.set_inlets(pores=self.net.pores('top'))
        assert sp.sum(self.alg['pore.inlets']) == 25

        self.alg.set_inlets(pores=self.net.pores('bottom'))
        assert sp.sum(self.alg['pore.inlets']) == 50

        self.alg.set_inlets(pores=self.net.pores('top'), overwrite=True)
        assert sp.sum(self.alg['pore.inlets']) == 25

        self.alg.set_inlets(pores=[], overwrite=True)
        assert sp.sum(self.alg['pore.inlets']) == 0

    def test_set_inlets_conflicting_with_outlets(self):
        self.alg = op.algorithms.Porosimetry(network=self.net)
        self.alg.setup(phase=self.water)
        self.alg['pore.outlets'][self.net.pores('top')] = True
        with pytest.raises(Exception):
            self.alg.set_inlets(pores=self.net.pores('top'))

    def test_set_outlets_conflicting_with_inlets(self):
        self.alg = op.algorithms.Porosimetry(network=self.net)
        self.alg.setup(phase=self.water)
        self.alg['pore.inlets'][self.net.pores('top')] = True
        with pytest.raises(Exception):
            self.alg.set_outlets(pores=self.net.pores('top'))

    def test_set_outlets_without_trapping(self):
        self.alg = op.algorithms.Porosimetry(network=self.net)
        self.alg.setup(phase=self.water)
        self.alg.set_inlets(pores=self.net.pores('top'))
        with pytest.raises(Exception):
            self.alg.set_outlets(pores=self.net.pores('top'))

    def test_set_outlets_overwrite(self):
        self.alg = op.algorithms.Porosimetry(network=self.net)
        self.alg.setup(phase=self.water)

        self.alg.set_outlets(pores=self.net.pores('top'))
        assert sp.sum(self.alg['pore.outlets']) == 25

        self.alg.set_outlets(pores=self.net.pores('bottom'))
        assert sp.sum(self.alg['pore.outlets']) == 50

        self.alg.set_outlets(pores=self.net.pores('top'), overwrite=True)
        assert sp.sum(self.alg['pore.outlets']) == 25

        self.alg.set_outlets(pores=[], overwrite=True)
        assert sp.sum(self.alg['pore.outlets']) == 0

    def test_set_residual_modes(self):
        self.alg = op.algorithms.Porosimetry(network=self.net)
        self.alg.setup(phase=self.water)

        Ps = sp.random.randint(0, self.net.Np, 10)
        Ts = self.net.find_neighbor_pores(pores=Ps)
        self.alg.set_residual(pores=Ps, throats=Ts)
        assert sp.sum(self.alg['pore.residual']) == sp.size(sp.unique(Ps))
        assert sp.sum(self.alg['throat.residual']) == sp.size(sp.unique(Ts))

        Ps = sp.random.randint(0, self.net.Np, 10)
        Ts = self.net.find_neighbor_pores(pores=Ps)
        self.alg.set_residual(pores=Ps, throats=Ts)
        assert sp.sum(self.alg['pore.residual']) > sp.size(sp.unique(Ps))
        assert sp.sum(self.alg['throat.residual']) > sp.size(sp.unique(Ts))

        Ps = sp.random.randint(0, self.net.Np, 10)
        Ts = self.net.find_neighbor_pores(pores=Ps)
        self.alg.set_residual(pores=Ps, throats=Ts, overwrite=True)
        assert sp.sum(self.alg['pore.residual']) == sp.size(sp.unique(Ps))
        assert sp.sum(self.alg['throat.residual']) == sp.size(sp.unique(Ts))

        self.alg.set_residual(pores=[], throats=[], overwrite=True)
        assert sp.sum(self.alg['pore.residual']) == 0

        self.alg.set_residual(pores=Ps, throats=Ts)
        self.alg.set_residual(overwrite=True)
        assert sp.sum(self.alg['pore.residual']) == 0

    def test_run_npts(self):
        self.alg = op.algorithms.Porosimetry(network=self.net)
        self.alg.setup(phase=self.water)
        Ps = sp.random.randint(0, self.net.Np, 10)
        self.alg.set_inlets(pores=Ps)
        self.alg.run(points=20)

    def test_run_inv_pressures(self):
        self.alg = op.algorithms.Porosimetry(network=self.net)
        self.alg.setup(phase=self.water)
        Ps = sp.random.randint(0, self.net.Np, 10)
        self.alg.set_inlets(pores=Ps)
        self.alg.run(points=range(0, 20000, 1000))

    def test_run_no_inlets(self):
        self.alg = op.algorithms.Porosimetry(network=self.net)
        self.alg.setup(phase=self.water)
        with pytest.raises(Exception):
            self.alg.run()

    def test_run_w_residual_pores_and_throats(self):
        self.alg = op.algorithms.Porosimetry(network=self.net)
        self.alg.setup(phase=self.water)
        self.alg.set_inlets(pores=self.net.pores('top'))
        self.alg.set_residual(pores=self.net.pores('bottom'))
        self.alg.run()
        data = self.alg.get_intrusion_data()
        assert hasattr(data, 'Pcap')
        assert hasattr(data, 'Snwp')


if __name__ == '__main__':

    t = PorosimetryTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
