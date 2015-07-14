import OpenPNM
import scipy as sp


class CapillaryPressureTest:

    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        self.geo = OpenPNM.Geometry.GenericGeometry(network=self.net,
                                                    pores=self.net.Ps,
                                                    throats=self.net.Ts)
        self.geo['throat.diameter'] = 1
        self.water = OpenPNM.Phases.GenericPhase(network=self.net)
        self.phys = OpenPNM.Physics.GenericPhysics(network=self.net,
                                                   geometry=self.geo,
                                                   phase=self.water)

    def test_washburn_pore_values(self):
        self.water['pore.surface_tension'] = 0.072
        self.water['pore.contact_angle'] = 120
        f = OpenPNM.Physics.models.capillary_pressure.washburn
        self.phys.models.add(propname='throat.capillary_pressure',
                             model=f,
                             surface_tension='pore.surface_tension',
                             contact_angle='pore.contact_angle',
                             throat_diameter='throat.diameter')
        a = 0.14399999999999993
        assert sp.allclose(self.water['throat.capillary_pressure'][0], a)
        self.phys.models.remove('throat.capillary_pressure')
        del self.water['pore.surface_tension']
        del self.water['pore.contact_angle']

    def test_washburn_throat_values(self):
        self.water['throat.surface_tension'] = 0.072
        self.water['throat.contact_angle'] = 120
        f = OpenPNM.Physics.models.capillary_pressure.washburn
        self.phys.models.add(propname='throat.capillary_pressure',
                             model=f,
                             surface_tension='throat.surface_tension',
                             contact_angle='throat.contact_angle',
                             throat_diameter='throat.diameter')
        a = 0.14399999999999993
        assert sp.allclose(self.water['throat.capillary_pressure'][0], a)
        self.phys.models.remove('throat.capillary_pressure')
        del self.water['throat.surface_tension']
        del self.water['throat.contact_angle']

    def test_purcell_pore_values(self):
        self.water['pore.surface_tension'] = 0.072
        self.water['pore.contact_angle'] = 120
        f = OpenPNM.Physics.models.capillary_pressure.purcell
        self.phys.models.add(propname='throat.capillary_pressure',
                             model=f,
                             r_toroid=0.1,
                             surface_tension='pore.surface_tension',
                             contact_angle='pore.contact_angle',
                             throat_diameter='throat.diameter')
        a = 0.26206427646507374
        assert sp.allclose(self.water['throat.capillary_pressure'][0], a)
        self.phys.models.remove('throat.capillary_pressure')
        del self.water['pore.surface_tension']
        del self.water['pore.contact_angle']

    def test_purcell_throat_values(self):
        self.water['throat.surface_tension'] = 0.072
        self.water['throat.contact_angle'] = 120
        f = OpenPNM.Physics.models.capillary_pressure.purcell
        self.phys.models.add(propname='throat.capillary_pressure',
                             model=f,
                             r_toroid=0.1,
                             surface_tension='throat.surface_tension',
                             contact_angle='throat.contact_angle',
                             throat_diameter='throat.diameter')
        a = 0.26206427646507374
        assert sp.allclose(self.water['throat.capillary_pressure'][0], a)
        self.phys.models.remove('throat.capillary_pressure')
        del self.water['throat.surface_tension']
        del self.water['throat.contact_angle']

    def test_static_pressure(self):
        self.water['pore.density'] = 1000
        self.water['pore.occupancy'] = True
        f = OpenPNM.Physics.models.capillary_pressure.static_pressure
        self.phys.models.add(propname='pore.static_pressure',
                             model=f,
                             gravity=[0, 0, 9.81],
                             pore_density='pore.density',
                             pore_occupancy='pore.occupancy')
        a = [0., 9810., 19620.]
        assert sp.all(sp.unique(self.water['pore.static_pressure']) == a)
