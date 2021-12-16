import numpy as np
import scipy as sp
import openpnm as op
from numpy.testing import assert_approx_equal


class VaporPressureTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.phase = op.phase.GenericPhase(network=self.net)
        self.phase['pore.temperature'] = 300*np.ones(self.phase.Np,)
        self.phase['pore.salinity'] = np.zeros((self.phase.Np,))

    def test_antoine(self):
        f = op.models.phase.vapor_pressure.antoine
        self.phase.add_model(propname='pore.test',
                             model=f,
                             temperature='pore.temperature',
                             A=8.088,
                             B=1750.71,
                             C=236.191)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.test'].mean(),
                            3607.8508)

    def test_water(self):
        f = op.models.phase.vapor_pressure.water
        self.phase.add_model(propname='pore.test',
                             model=f,
                             temperature='pore.temperature',
                             salinity='pore.salinity')
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.test'].mean(),
                            3536.0130)

    def test_water_no_salinity(self):
        f = op.models.phase.vapor_pressure.water
        del self.phase['pore.salinity']
        self.phase.add_model(propname='pore.test',
                             model=f,
                             temperature='pore.temperature',
                             salinity='pore.salinity')
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.test'].mean(),
                            3536.0130)
        self.phase['pore.salinity'] = np.zeros((self.phase.Np,))


if __name__ == '__main__':

    t = VaporPressureTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
