import pytest
import numpy as np
import openpnm as op
mgr = op.Workspace()


class OrdinaryPercolationTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5], spacing=0.0005)
        self.geo = op.geometry.SpheresAndCylinders(network=self.net,
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
        self.alg = op.algorithms.OrdinaryPercolation(network=self.net,
                                                     phase=self.water)

        self.alg.set_inlets(pores=self.net.pores('top'))
        assert np.sum(self.alg['pore.inlets']) == 25

        self.alg.set_inlets(pores=self.net.pores('bottom'))
        assert np.sum(self.alg['pore.inlets']) == 50

        self.alg.set_inlets(pores=self.net.pores('top'), overwrite=True)
        assert np.sum(self.alg['pore.inlets']) == 25

        self.alg.set_inlets(pores=[], overwrite=True)
        assert np.sum(self.alg['pore.inlets']) == 0

    def test_set_inlets_conflicting_with_outlets(self):
        self.alg = op.algorithms.OrdinaryPercolation(network=self.net,
                                                     phase=self.water)
        self.alg['pore.outlets'][self.net.pores('top')] = True
        with pytest.raises(Exception):
            self.alg.set_inlets(pores=self.net.pores('top'))

    def test_set_outlets_conflicting_with_inlets(self):
        self.alg = op.algorithms.OrdinaryPercolation(network=self.net,
                                                     phase=self.water)
        self.alg['pore.inlets'][self.net.pores('top')] = True
        with pytest.raises(Exception):
            self.alg.set_outlets(pores=self.net.pores('top'))

    def test_set_outlets_without_trapping(self):
        self.alg = op.algorithms.OrdinaryPercolation(network=self.net,
                                                     phase=self.water)
        self.alg.set_inlets(pores=self.net.pores('top'))
        with pytest.raises(Exception):
            self.alg.set_outlets(pores=self.net.pores('top'))

    def test_set_outlets_overwrite(self):
        self.alg = op.algorithms.OrdinaryPercolation(network=self.net,
                                                     phase=self.water)

        self.alg.set_outlets(pores=self.net.pores('top'))
        assert np.sum(self.alg['pore.outlets']) == 25

        self.alg.set_outlets(pores=self.net.pores('bottom'))
        assert np.sum(self.alg['pore.outlets']) == 50

        self.alg.set_outlets(pores=self.net.pores('top'), overwrite=True)
        assert np.sum(self.alg['pore.outlets']) == 25

        self.alg.set_outlets(pores=[], overwrite=True)
        assert np.sum(self.alg['pore.outlets']) == 0

    def test_set_residual_modes(self):
        self.alg = op.algorithms.OrdinaryPercolation(network=self.net,
                                                     phase=self.water)

        Ps = np.random.randint(0, self.net.Np, 10)
        Ts = self.net.find_neighbor_pores(pores=Ps)
        self.alg.set_residual(pores=Ps, throats=Ts)
        assert np.sum(self.alg['pore.residual']) == np.size(np.unique(Ps))
        assert np.sum(self.alg['throat.residual']) == np.size(np.unique(Ts))

        Ps = np.random.randint(0, self.net.Np, 10)
        Ts = self.net.find_neighbor_pores(pores=Ps)
        self.alg.set_residual(pores=Ps, throats=Ts)
        assert np.sum(self.alg['pore.residual']) > np.size(np.unique(Ps))
        assert np.sum(self.alg['throat.residual']) > np.size(np.unique(Ts))

        Ps = np.random.randint(0, self.net.Np, 10)
        Ts = self.net.find_neighbor_pores(pores=Ps)
        self.alg.set_residual(pores=Ps, throats=Ts, overwrite=True)
        assert np.sum(self.alg['pore.residual']) == np.size(np.unique(Ps))
        assert np.sum(self.alg['throat.residual']) == np.size(np.unique(Ts))

        self.alg.set_residual(pores=[], throats=[], overwrite=True)
        assert np.sum(self.alg['pore.residual']) == 0

        self.alg.set_residual(pores=Ps, throats=Ts)
        self.alg.set_residual(overwrite=True)
        assert np.sum(self.alg['pore.residual']) == 0

    def test_run_npts(self):
        self.alg = op.algorithms.OrdinaryPercolation(network=self.net,
                                                     phase=self.water)
        Ps = np.random.randint(0, self.net.Np, 10)
        self.alg.set_inlets(pores=Ps)
        self.alg.run(points=20)

    def test_run_inv_pressures(self):
        self.alg = op.algorithms.OrdinaryPercolation(network=self.net,
                                                     phase=self.water)
        Ps = np.random.randint(0, self.net.Np, 10)
        self.alg.set_inlets(pores=Ps)
        self.alg.run(points=range(0, 20000, 1000))

    def test_run_no_inlets(self):
        self.alg = op.algorithms.OrdinaryPercolation(network=self.net,
                                                     phase=self.water)
        with pytest.raises(Exception):
            self.alg.run()

    def test_run_w_residual_pores_and_throats(self):
        self.alg = op.algorithms.OrdinaryPercolation(network=self.net,
                                                     phase=self.water)
        self.alg.set_inlets(pores=self.net.pores('top'))
        self.alg.set_residual(pores=self.net.pores('bottom'))
        self.alg.run()
        data = self.alg.results(Pc=20000)
        assert sum(data['pore.occupancy']) > 0
        assert sum(data['throat.occupancy']) > 0

    def test_is_percolating(self):
        self.alg = op.algorithms.OrdinaryPercolation(network=self.net,
                                                     phase=self.water)
        self.alg.settings['access_limited'] = True
        self.alg.set_inlets(pores=self.net.pores('top'))
        self.alg.set_outlets(pores=self.net.pores('bottom'))
        self.alg.run()
        assert not self.alg.is_percolating(0)
        assert self.alg.is_percolating(1e5)

    def test_entry_vs_invasion_pressure(self):
        self.alg = op.algorithms.OrdinaryPercolation(network=self.net,
                                                     phase=self.water)
        self.alg.settings._update({'mode': 'bond',
                                   'access_limited': True})
        self.alg.set_inlets(pores=self.net.pores('top'))
        self.alg.set_outlets(pores=self.net.pores('bottom'))
        self.alg.run()
        Tinv = self.alg['throat.invasion_pressure']
        Tent = self.water['throat.entry_pressure']
        assert np.all(Tent <= Tinv)


if __name__ == '__main__':

    t = OrdinaryPercolationTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
