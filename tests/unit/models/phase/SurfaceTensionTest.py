import openpnm as op
from numpy.testing import assert_approx_equal, assert_allclose
import chemicals


class SurfaceTensionTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.phase = op.phase.Species(network=self.net, species='water')
        self.phase['pore.salinity'] = 0  # g/kg

    def test_water(self):
        f = op.models.phase.surface_tension.water_correlation
        self.phase.add_model(propname='pore.surface_tension',
                             model=f)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.surface_tension'].mean(),
                            0.07199533)

    def test_brock_bird(self):
        f = op.models.phase.surface_tension.liquid_pure_brock_bird
        self.phase.add_model(propname='pore.surface_tension',
                             model=f)
        self.phase.regenerate_models()
        # assert_approx_equal(self.phase['pore.surface_tension'].mean(),
        #                     0.07820759)

    def test_generic_chemicals_for_liquid_mixtures(self):
        h2o = op.phase.Species(network=self.net, species='water')
        h2o.add_model(propname='pore.surface_tension',
                      model=op.models.phase.surface_tension.liquid_pure_brock_bird)
        etoh = op.phase.Species(network=self.net, species='ethanol')
        etoh.add_model(propname='pore.surface_tension',
                       model=op.models.phase.surface_tension.liquid_pure_brock_bird)
        vodka = op.phase.LiquidMixture(network=self.net, components=[h2o, etoh])
        vodka.x(h2o.name, 0.5)
        vodka.x(etoh.name, 0.5)
        mods = [
            chemicals.interface.Winterfeld_Scriven_Davis
        ]
        vals = []
        for f in mods:
            vals.append(op.models.phase.chemicals_wrapper(target=vodka, f=f).mean())
        assert_allclose(vals, 2.898e-1, rtol=.8)

    def test_generic_chemicals_for_pure_liq(self):
        mods = [
            # chemicals.interface.REFPROP_sigma,  # Needs sigma0
            # chemicals.interface.Somayajulu,  # Needs A
            # chemicals.interface.Jasper,  # Needs a
            chemicals.interface.Brock_Bird,
            # chemicals.interface.Sastri_Rao,  # Numba version not working
            # chemicals.interface.Pitzer,  # Model missing
            chemicals.interface.Zuo_Stenby,
            chemicals.interface.Miqueu,
            # chemicals.interface.Aleem,  # Needs rhol
            chemicals.interface.sigma_IAPWS,
        ]
        h2o = op.phase.Species(network=self.net, species='water')
        vals = []
        for f in mods:
            vals.append(op.models.phase.chemicals_wrapper(target=h2o, f=f).mean())
        assert_allclose(vals, 2.898e-1, rtol=.8)


if __name__ == '__main__':

    t = SurfaceTensionTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
