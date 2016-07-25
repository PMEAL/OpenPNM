import OpenPNM
import scipy as sp


class VaporPressureTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        self.phase = OpenPNM.Phases.GenericPhase(network=self.net)
        self.phase['pore.temperature'] = 300*sp.ones(self.phase.Np,)
        self.phase['pore.salinity'] = sp.zeros((self.phase.Np,))

    def test_antoine(self):
        f = OpenPNM.Phases.models.vapor_pressure.antoine
        self.phase.models.add(propname='pore.test',
                              model=f,
                              pore_temperature='pore.temperature',
                              A=8.07131,
                              B=1730.63,
                              C=233.426)
        assert sp.allclose(self.phase['pore.test'], 3523.72641773)

    def test_water(self):
        f = OpenPNM.Phases.models.vapor_pressure.water
        self.phase.models.add(propname='pore.test',
                              model=f,
                              pore_temperature='pore.temperature',
                              pore_salinity='pore.salinity')
        assert sp.allclose(self.phase['pore.test'], 3536.01)

    def test_water_no_salinity(self):
        f = OpenPNM.Phases.models.vapor_pressure.water
        del self.phase['pore.salinity']
        self.phase.models.add(propname='pore.test',
                              model=f,
                              pore_temperature='pore.temperature',
                              pore_salinity='pore.salinity')
        assert sp.allclose(self.phase['pore.test'], 3536.01)
        self.phase['pore.salinity'] = sp.zeros((self.phase.Np,))
