import OpenPNM
import scipy as sp


class SurfaceTensionTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        self.phase = OpenPNM.Phases.GenericPhase(network=self.net)
        self.phase['pore.temperature'] = 298.0  # K
        self.phase['pore.molecular_weight'] = 0.018  # kg/mol
        self.phase['pore.critical_temperature'] = 647.15  # K
        self.phase['pore.critical_pressure'] = 3771000.0  # Pa
        self.phase['pore.salinity'] = 0  # g/kg
        self.phase['pore.molar_density'] = 55.5  # mol/m3

    def test_water(self):
        f = OpenPNM.Phases.models.surface_tension.water
        self.phase.models.add(propname='pore.surface_tension',
                              model=f)
        assert sp.allclose(self.phase['pore.surface_tension'], 0.07199533)

    def test_eotvos(self):
        f = OpenPNM.Phases.models.surface_tension.eotvos
        self.phase.models.add(propname='pore.surface_tension',
                              model=f,
                              k=0.000014)
        assert sp.allclose(self.phase['pore.surface_tension'], 0.07112169)

    def test_guggenheim_katayama(self):
        f = OpenPNM.Phases.models.surface_tension.guggenheim_katayama
        self.phase.models.add(propname='pore.surface_tension',
                              model=f,
                              K2=0.0000014,
                              n=0.1)
        assert sp.allclose(self.phase['pore.surface_tension'], 0.27582571)

    def test_brock_bird_scaling(self):
        f = OpenPNM.Phases.models.surface_tension.brock_bird_scaling
        self.phase.models.add(propname='pore.surface_tension',
                              model=f,
                              sigma_o=0.0608,
                              To=363)
        assert sp.allclose(self.phase['pore.surface_tension'], 0.07820761)
