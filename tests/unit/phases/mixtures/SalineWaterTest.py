import openpnm as op
from openpnm.phase import mixtures


class SalineWaterTest:
    def setup_class(self):
        ws = op.Workspace()
        ws.clear()
        self.net = op.network.Cubic(shape=[10, 10, 10])

    def test_init(self):
        self.sw = mixtures.SalineWater(network=self.net)
        assert isinstance(self.sw, mixtures.GenericMixture)


if __name__ == '__main__':

    t = SalineWaterTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
