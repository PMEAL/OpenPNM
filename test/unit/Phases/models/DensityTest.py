import OpenPNM
import scipy as sp


class DensityTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        self.phase = OpenPNM.Phases.GenericPhase(network=self.net)
        self.phase['pore.temperature'] = 298.0  # K
        self.phase['pore.pressure'] = 101325.0  # Pa
        self.phase['pore.molecular_weight'] = 0.018  # kg/mol
        self.phase['pore.molar_density'] = 55539.0  # mol/m3
        self.phase['pore.salinity'] = 0.0  # ppt

    def test_standard(self):
        # Liquid water
        self.phase.models.add(propname='pore.density',
                              model=OpenPNM.Phases.models.density.standard)
        assert sp.allclose(self.phase['pore.density'], 999.702)

    def test_ideal_gas(self):
        # Water vapor
        self.phase.models.add(propname='pore.density',
                              model=OpenPNM.Phases.models.density.ideal_gas)
        assert sp.allclose(self.phase['pore.density'], 0.73610248)

    def test_water(self):
        # Liquid water
        self.phase.models.add(propname='pore.density',
                              model=OpenPNM.Phases.models.density.water)
        assert sp.allclose(self.phase['pore.density'], 996.9522)

    def teardown_class(self):
        del(self.phase)
        del(self.net)
