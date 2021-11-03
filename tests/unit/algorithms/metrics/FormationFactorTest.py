import numpy as np
import openpnm as op
mgr = op.Workspace()


class FormationFactorTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[15, 15, 15], spacing=0.0005)
        self.geo = op.geometry.SpheresAndCylinders(network=self.net,
                                                   pores=self.net.Ps,
                                                   throats=self.net.Ts)

    def test_run(self):
        FF = op.algorithms.metrics.FormationFactor(network=self.net)
        assert len(FF.results) == 0
        FF.run()
        assert len(FF.results) == 3

    def test_given_area(self):
        FF = op.algorithms.metrics.FormationFactor(network=self.net)
        FF.run()
        val_1 = FF.results['x']
        FF.set_area(direction='x', area=(15*0.0005)**2)
        FF.run()
        val_2 = FF.results['x']
        assert val_1 != val_2

    def test_given_length(self):
        FF = op.algorithms.metrics.FormationFactor(network=self.net)
        FF.run()
        val_1 = FF.results['y']
        FF.set_length(direction='y', length=15*0.0005)
        FF.run()
        val_2 = FF.results['y']
        assert val_1 != val_2

    def test_setting_inlets(self):
        FF = op.algorithms.metrics.FormationFactor(network=self.net)
        Ps = self.net.pores('left')
        self.net.set_label(pores=Ps, label='blah')
        FF.set_inlets(direction='x', label='blah')
        FF.run()
        val_1 = FF.results['x']
        FF.set_inlets(direction='x', label='left')
        FF.run()
        val_2 = FF.results['x']
        np.testing.assert_allclose(val_1, val_2)


if __name__ == '__main__':

    t = FormationFactorTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('Running test: '+item)
            t.__getattribute__(item)()
