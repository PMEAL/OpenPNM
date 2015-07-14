import OpenPNM
import scipy as sp


class ElectricalConductanceTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[4, 4, 4])
        self.geo = OpenPNM.Geometry.GenericGeometry(network=self.net,
                                                    pores=self.net.Ps,
                                                    throats=self.net.Ts)
        self.phase = OpenPNM.Phases.GenericPhase(network=self.net)
        self.phys = OpenPNM.Physics.GenericPhysics(network=self.net,
                                                   phase=self.phase,
                                                   geometry=self.geo)

    def test_electrical_conductance(self):
        self.phase['pore.conductivity'] = 1
        self.geo['pore.area'] = 1
        self.geo['pore.diameter'] = 1
        self.geo['throat.area'] = 1
        self.geo['throat.length'] = 0.0001
        f = OpenPNM.Physics.models.electrical_conductance.series_resistors
        self.phys.models.add(propname='throat.electrical_conductance',
                             pore_conductivity='pore.conductivity',
                             model=f)
        a = 0.99990001
        assert sp.allclose(self.phys['throat.electrical_conductance'][0], a)
