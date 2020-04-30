import openpnm as op
import openpnm.models.geometry.throat_length as mods
import numpy as np
from numpy.testing import assert_approx_equal


class ThroatStraightTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[1, 2, 1], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.base = np.array([0.5, 0.5, 0.5])
        self.geo['throat.endpoints.head'] = np.array([[0, 0.2, 0]]) + self.base
        self.geo['throat.endpoints.tail'] = np.array([[0, 0.7, 0]]) + self.base

    def test_piecewise(self):
        self.geo.add_model(propname='throat.length',
                           model=mods.piecewise,
                           regen_mode='normal')
        actual = self.geo['throat.length'][0]
        assert_approx_equal(actual, desired=0.5)

    def test_piecewise_with_centroid(self):
        self.geo['throat.centroid'] = np.array([[0, 0.5, 0.5]]) + self.base
        self.geo.add_model(propname='throat.length',
                           model=mods.piecewise,
                           regen_mode='normal')
        actual = self.geo['throat.length'][0]
        assert_approx_equal(actual, desired=1.1216117)

    def test_conduit_lengths_wo_centroid(self):
        del self.geo["throat.centroid"]
        self.geo['throat.length'] = 0.5
        self.geo.add_model(propname='throat.conduit_lengths',
                           model=mods.conduit_lengths,
                           regen_mode='normal')
        L1 = self.geo['throat.conduit_lengths.pore1'][0]
        L2 = self.geo['throat.conduit_lengths.pore2'][0]
        Lt = self.geo['throat.conduit_lengths.throat'][0]
        assert_approx_equal(actual=L1, desired=0.2)
        assert_approx_equal(actual=L2, desired=0.3)
        assert_approx_equal(actual=Lt, desired=0.5)

    def test_conduit_lengths_w_centroid(self):
        # Delete "throat.length" key, otherwise, it won't get re-calculated
        del self.geo['throat.length']
        self.geo["throat.centroid"] = np.array([[0, 0.5, 0.5]]) + self.base
        self.geo.add_model(propname='throat.conduit_lengths',
                           model=mods.conduit_lengths,
                           regen_mode='normal')
        L1 = self.geo['throat.conduit_lengths.pore1'][0]
        L2 = self.geo['throat.conduit_lengths.pore2'][0]
        Lt = self.geo['throat.conduit_lengths.throat'][0]
        assert_approx_equal(actual=L1, desired=0.2)
        assert_approx_equal(actual=L2, desired=0.3)
        assert_approx_equal(actual=Lt, desired=1.12161167)


if __name__ == '__main__':

    t = ThroatStraightTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
