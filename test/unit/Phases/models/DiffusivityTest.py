import OpenPNM
import scipy as sp


class DiffusivityTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        self.phase = OpenPNM.Phases.GenericPhase(network=self.net)
        self.phase['pore.temperature'] = 298.0  # K
        self.phase['pore.pressure'] = 101325  # Pa
        self.phase['pore.viscosity'] = 1.75e-5  # Pa.s

    def test_fuller(self):
        f = OpenPNM.Phases.models.diffusivity.fuller
        self.phase.models.add(propname='pore.diffusivity',
                              model=f,
                              MA=0.032,
                              MB=0.028,
                              vA=16.6,
                              vB=17.9)
        assert sp.allclose(self.phase['pore.diffusivity'], 2.06754784e-05)

    def test_fuller_scaling(self):
        f = OpenPNM.Phases.models.diffusivity.fuller_scaling
        self.phase.models.add(propname='pore.diffusivity',
                              model=f,
                              DABo=1.79712526e-05,
                              Po=100000,
                              To=273)
        assert sp.allclose(self.phase['pore.diffusivity'], 2.06754784e-05)

    def test_tyn_calus(self):
        f = OpenPNM.Phases.models.diffusivity.tyn_calus
        self.phase.models.add(propname='pore.diffusivity',
                              model=f,
                              VA=16.5,
                              VB=17.9,
                              sigma_A=1,
                              sigma_B=1)
        assert sp.allclose(self.phase['pore.diffusivity'], 9.84851806e-05)

    def test_tyn_calus_scaling(self):
        f = OpenPNM.Phases.models.diffusivity.tyn_calus_scaling
        self.phase.models.add(propname='pore.diffusivity',
                              model=f,
                              DABo=5.26300839e-05,
                              mu_o=3e-5,
                              To=273)
        assert sp.allclose(self.phase['pore.diffusivity'], 9.84851806e-05)
