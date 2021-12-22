import openpnm as op
import scipy as sp
from numpy.testing import assert_approx_equal


class SurfaceTensionTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.phase = op.phase.GenericPhase(network=self.net)
        self.phase['pore.temperature'] = 298.0  # K
        self.phase['pore.molecular_weight'] = 0.018  # kg/mol
        self.phase['pore.critical_temperature'] = 647.15  # K
        self.phase['pore.critical_pressure'] = 3771000.0  # Pa
        self.phase['pore.salinity'] = 0  # g/kg
        self.phase['pore.molar_density'] = 55.5  # mol/m3

    def test_water(self):
        f = op.models.phase.surface_tension.water
        self.phase.add_model(propname='pore.surface_tension',
                             model=f)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.surface_tension'].mean(),
                            0.07199533)

    def test_eotvos(self):
        f = op.models.phase.surface_tension.eotvos
        self.phase.add_model(propname='pore.surface_tension',
                             model=f,
                             k=0.000014)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.surface_tension'].mean(),
                            0.07112169)

    def test_guggenheim_katayama(self):
        f = op.models.phase.surface_tension.guggenheim_katayama
        self.phase.add_model(propname='pore.surface_tension',
                             model=f,
                             K2=0.0000014,
                             n=0.1)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.surface_tension'].mean(),
                            0.27582571)

    def test_brock_bird_scaling(self):
        f = op.models.phase.surface_tension.brock_bird_scaling
        self.phase.add_model(propname='pore.surface_tension',
                             model=f,
                             sigma_o=0.0608,
                             To=363)
        self.phase.regenerate_models()
        assert_approx_equal(self.phase['pore.surface_tension'].mean(),
                            0.07820759)


if __name__ == '__main__':

    t = SurfaceTensionTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
