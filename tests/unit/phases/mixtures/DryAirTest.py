import openpnm as op
from openpnm.phase import mixtures


class DryAirTest:
    def setup_class(self):
        ws = op.Workspace()
        ws.clear()
        self.net = op.network.Cubic(shape=[10, 10, 10])

    def test_init(self):
        da = mixtures.DryAir(network=self.net)
        assert isinstance(da, mixtures.GenericMixture)


if __name__ == '__main__':

    t = DryAirTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
