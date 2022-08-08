import numpy as np
import openpnm as op
from numpy.testing import assert_approx_equal, assert_array_almost_equal, assert_allclose
import chemicals


class VaporPressureTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.phase = op.phase.Species(network=self.net, species='h2o')
        self.phase['pore.temperature'] = 300*np.ones(self.phase.Np,)
        self.phase['pore.salinity'] = np.zeros((self.phase.Np,))

    def test_antoine(self):
        f = op.models.phase.vapor_pressure.liquid_pure_antoine
        self.phase.add_model(propname='pore.test',
                             model=f,
                             T='pore.temperature')
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.test'].mean(), 3546.9831)

    def test_water(self):
        f = op.models.phase.vapor_pressure.water_correlation
        self.phase.add_model(propname='pore.test',
                             model=f,
                             T='pore.temperature',
                             salinity='pore.salinity')
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.test'].mean(), 3536.0130)

    def test_water_no_salinity(self):
        f = op.models.phase.vapor_pressure.water_correlation
        del self.phase['pore.salinity']
        self.phase.add_model(propname='pore.test',
                             model=f,
                             T='pore.temperature',
                             salinity='pore.salinity')
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.test'].mean(), 3536.0130)
        self.phase['pore.salinity'] = np.zeros((self.phase.Np,))

    def test_generic_chemicals_for_pure_liquid(self):
        mods = [
            # chemicals.vapor_pressure.Wagner_original,  # Needs constants
            # chemicals.vapor_pressure.Wagner,  # Needs constants
            # chemicals.vapor_pressure.TRC_Antoine_extended,  # Needs constants
            # chemicals.vapor_pressure.Antoine,  # Needs constants
            chemicals.vapor_pressure.boiling_critical_relation,
            chemicals.vapor_pressure.Lee_Kesler,
            chemicals.vapor_pressure.Ambrose_Walton,
            chemicals.vapor_pressure.Sanjari,
            chemicals.vapor_pressure.Edalat,
            chemicals.iapws.iapws95_Psat,
        ]
        h2o = op.phase.Species(network=self.net, species='water')
        vals = []
        for f in mods:
            vals.append(op.models.phase.chemicals_wrapper(phase=h2o, f=f).mean())
        assert_allclose(vals, 2.762e3, rtol=.5)


if __name__ == '__main__':

    t = VaporPressureTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
