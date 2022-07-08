import openpnm as op
import openpnm.models.phase as pm
from numpy.testing import assert_allclose


class SpeciesTest:

    def setup_class(self):
        pass

    def test_liquid_by_name(self):
        pn = op.network.Demo()
        h2o = op.phase.Species(network=pn, species='water')
        h2o.add_model(propname='pore.density',
                      model=pm.density.liquid_density)
        h2o.add_model(propname='pore.vapor_pressure',
                      model=pm.vapor_pressure.vapor_pressure)
        h2o.add_model(propname='pore.thermal_conductivity',
                      model=pm.thermal_conductivity.liquid_thermal_conductivity)
        assert_allclose(h2o['pore.thermal_conductivity'].mean(), 0.658328)
        assert_allclose(h2o['pore.density'].mean(), 992.34552)
        assert_allclose(h2o['pore.vapor_pressure'], 3150.40831432)


if __name__ == '__main__':

    t = SpeciesTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
