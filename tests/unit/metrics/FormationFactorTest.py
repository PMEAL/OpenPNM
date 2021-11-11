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
        FF = op.metrics.FormationFactor(network=self.net)
        F = FF.run()
        np.testing.assert_allclose(F, 22.37677317)

    def test_given_area(self):
        FF = op.metrics.FormationFactor(network=self.net)
        val_1 = FF.run()
        FF.settings.update({'area': (15*0.0005)**2})
        val_2 = FF.run()
        assert val_1 != val_2

    def test_given_length(self):
        FF = op.metrics.FormationFactor(network=self.net)
        val_1 = FF.run()
        FF.settings.update({'length': 15*0.0005})
        val_2 = FF.run()
        assert val_1 != val_2

    def test_setting_inlets(self):
        FF = op.metrics.FormationFactor(network=self.net)
        FF.settings.update({'inlet': 'top', 'outlet': 'bottom'})
        val_1 = FF.run()
        FF.settings.update({'inlet': 'front', 'outlet': 'back'})
        val_2 = FF.run()
        assert val_1 != val_2


if __name__ == '__main__':

    t = FormationFactorTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('Running test: '+item)
            t.__getattribute__(item)()
