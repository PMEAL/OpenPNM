import openpnm as op
import scipy as sp
import pytest


class ModelsTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)

    def teardown_class(self):
        ws = op.core.Workspace()
        ws.clear()

    def test_models_dict_print(self):
        s = self.geo.models.__str__().split('\n')
        assert len(s) == 55
        assert s.count('―'*78) == 13

    def test_regenerate_models(self):
        a = len(self.geo.props())
        assert a == 11
        self.geo.clear(mode='props')
        a = len(self.geo.props())
        assert a == 0


if __name__ == '__main__':

    t = ModelsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
