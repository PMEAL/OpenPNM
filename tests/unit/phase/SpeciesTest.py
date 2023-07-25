from numpy.testing import assert_allclose
from thermo import Chemical

import openpnm as op
import openpnm.models.phase as pm


class SpeciesTest:

    def setup_class(self):
        pass

    def test_species_by_name(self):
        pn = op.network.Demo()
        h2o = op.phase.Species(network=pn, species='water')
        h2o.add_model(propname='pore.density',
                      model=pm.density.liquid_pure_COSTALD)
        h2o.add_model(propname='pore.vapor_pressure',
                      model=pm.vapor_pressure.liquid_pure_antoine)
        h2o.add_model(propname='pore.thermal_conductivity',
                      model=pm.thermal_conductivity.liquid_pure_gismr)
        assert_allclose(h2o['pore.thermal_conductivity'].mean(), 0.65836, rtol=1e-5)
        assert_allclose(h2o['pore.density'].mean(), 993.327626, rtol=1e-5)
        assert_allclose(h2o['pore.vapor_pressure'], 3150.40831432)

    def test_standard_liquid(self):
        pn = op.network.Demo()
        h2o = op.phase.StandardLiquid(network=pn, species='water')
        h2o.regenerate_models()
        a = Chemical('h2o')
        assert_allclose(h2o['pore.density'].mean(), a.rho, rtol=0.01)
        assert_allclose(h2o['pore.heat_capacity'].mean(), a.Cplm, rtol=0.5)
        assert_allclose(h2o['pore.thermal_conductivity'].mean(), a.kl, rtol=0.2)
        assert_allclose(h2o['pore.viscosity'].mean(), a.mul, rtol=0.5)
        assert_allclose(h2o['pore.vapor_pressure'].mean(), a.Psat, rtol=0.01)

    def test_standard_gas(self):
        pn = op.network.Demo()
        o2 = op.phase.StandardGas(network=pn, species='o2')
        o2.regenerate_models()
        a = Chemical('o2')
        assert_allclose(o2['pore.density'].mean(), a.rhog, rtol=0.01)
        assert_allclose(o2['pore.heat_capacity'].mean(), a.Cpgm, rtol=0.01)
        assert_allclose(o2['pore.thermal_conductivity'].mean(), a.kg, rtol=0.5)
        # assert_allclose(o2['pore.viscosity'].mean(), a.mug, rtol=0.5)

    def test_get_mixture(self):
        net = op.network.Demo()
        o2 = op.phase.StandardGas(network=net, species='o2', name='pure_O2')
        n2 = op.phase.StandardGas(network=net, species='n2', name='pure_N2')
        air = op.phase.GasMixture(network=net, components=[n2, o2])
        assert air is o2.mixture
        assert air is n2.mixture


if __name__ == '__main__':

    t = SpeciesTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f"Running test: {item}")
            t.__getattribute__(item)()
