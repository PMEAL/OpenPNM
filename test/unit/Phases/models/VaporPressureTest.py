import OpenPNM
import scipy as sp

class VaporPressureTest:
    def class_setup(self):
        self.net = OpenPNM.Network.Cubic(shape=[3,3,3])
        self.phase = OpenPNM.Phases.GenericPhase(network=self.net)
        self.phase['pore.temperature'] = 300*sp.ones(self.phase.Np,)
        self.phase['pore.salinity'] = sp.zeros((self.phase.Np,))

    def test_antoine(self):
        self.phase.models.add(propname='pore.test',
                              model=OpenPNM.Phases.models.vapor_pressure.antoine,
                              pore_temperature='pore.temperature',
                              A = 6.20963,
                              B = 2354.731,
                              C = 7.559)
        assert sp.allclose(self.phase['pore.test'], 0.0357632)

    def test_water(self):
        self.phase.models.add(propname='pore.test',
                              model=OpenPNM.Phases.models.vapor_pressure.water,
                              pore_temperature='pore.temperature',
                              pore_salinity='pore.salinity')
        assert sp.allclose(self.phase['pore.test'], 3536.01)