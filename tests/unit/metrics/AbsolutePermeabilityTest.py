import numpy as np
import openpnm as op
mgr = op.Workspace()


class FormationFactorTest:

    def setup_class(self):
        np.random.seed(5)
        self.net = op.network.Cubic(shape=[15, 15, 15], spacing=0.0005)
        self.geo = op.geometry.SpheresAndCylinders(network=self.net,
                                                   pores=self.net.Ps,
                                                   throats=self.net.Ts)

    def test_run(self):
        perm = op.metrics.AbsolutePermeability(network=self.net)
        K = perm.run()
        print(K)
        np.testing.assert_allclose(K, 1.71917212e-11)

    def test_given_area(self):
        perm = op.metrics.AbsolutePermeability(network=self.net)
        val_1 = perm.run()
        perm.settings._update({'area': (15*0.0005)**2})
        val_2 = perm.run()
        assert val_1 != val_2

    def test_given_length(self):
        perm = op.metrics.AbsolutePermeability(network=self.net)
        val_1 = perm.run()
        perm.settings._update({'length': 15*0.0005})
        val_2 = perm.run()
        assert val_1 != val_2

    def test_setting_inlets(self):
        perm = op.metrics.AbsolutePermeability(network=self.net)
        perm.settings._update({'inlet': 'top', 'outlet': 'bottom'})
        val_1 = perm.run()
        perm.settings._update({'inlet': 'front', 'outlet': 'back'})
        val_2 = perm.run()
        assert val_1 != val_2


if __name__ == '__main__':

    t = FormationFactorTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('Running test: '+item)
            t.__getattribute__(item)()
