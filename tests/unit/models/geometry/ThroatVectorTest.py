import openpnm as op
import openpnm.models.geometry as gm
import numpy as np


class ThroatVectorTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 1, 1], spacing=1)

    def test_pore_to_pore(self):
        self.net.add_model(propname='throat.vector',
                           model=gm.throat_vector.pore_to_pore)
        assert np.all(self.net['throat.vector'][:, 0] == 1.0)
        assert np.all(self.net['throat.vector'][:, 1] == 0.0)
        assert np.all(self.net['throat.vector'][:, 2] == 0.0)


if __name__ == '__main__':

    t = ThroatVectorTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
