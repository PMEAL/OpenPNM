import openpnm as op
from openpnm.utils import get_mixture_model_args
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

    def test_generic_chemicals_for_liquid_mixtures(self):
        pn = op.network.Demo()
        A = op.phase.Species(network=pn, species='chlorobenzene')
        A.add_model(propname='pore.surface_tension',
                    model=op.models.phase.surface_tension.liquid_pure_bb)
        A.add_model(propname='pore.density',
                    model=op.models.phase.density.liquid_pure_COSTALD)
        A['pore.molar_density'] = A['pore.density']/(A['param.molecular_weight']/1000)
        s_calc = A['pore.surface_tension'][0]
        s_ref = chemicals.interface.Brock_Bird(
            T=A['pore.temperature'][0],
            Tb=A['param.boiling_temperature'],
            Tc=A['param.critical_temperature'],
            Pc=A['param.critical_pressure'],
        )
        assert_allclose(s_ref, s_calc, rtol=1e-10, atol=0)

        B = op.phase.Species(network=pn, species='benzene')
        B.add_model(propname='pore.surface_tension',
                    model=op.models.phase.surface_tension.liquid_pure_bb)
        B.add_model(propname='pore.density',
                    model=op.models.phase.density.liquid_pure_COSTALD)
        B['pore.molar_density'] = B['pore.density']/(B['param.molecular_weight']/1000)
        s_calc = B['pore.surface_tension'][0]
        s_ref = chemicals.interface.Brock_Bird(
            T=B['pore.temperature'][0],
            Tb=B['param.boiling_temperature'],
            Tc=B['param.critical_temperature'],
            Pc=B['param.critical_pressure'],
        )
        assert_allclose(s_ref, s_calc, rtol=1e-10, atol=0)

        vodka = op.phase.LiquidMixture(network=pn, components=[A, B])
        vodka.x(A, 0.5)
        vodka.x(B, 0.5)
        vodka.add_model(propname='pore.surface_tension',
                        model=op.models.phase.surface_tension.liquid_mixture_wsd)
        s_calc = vodka['pore.surface_tension'][0]
        args = get_mixture_model_args(
            target=vodka,
            composition='xs',
            args={
                'sigmas': 'pore.surface_tension.*',
                'rhoms': 'pore.molar_density.*',
            })
        s_ref = chemicals.interface.Winterfeld_Scriven_Davis(**args)
        assert_allclose(s_ref, s_calc, rtol=1e-2, atol=0)

    def test_generic_chemicals_for_pure_liq(self):
        mods = [
            # chemicals.interface.REFPROP_sigma,  # Needs sigma0
            # chemicals.interface.Somayajulu,  # Needs A
            # chemicals.interface.Jasper,  # Needs a
            chemicals.interface.Brock_Bird,
            chemicals.interface.Sastri_Rao,  # Numba version not working
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
