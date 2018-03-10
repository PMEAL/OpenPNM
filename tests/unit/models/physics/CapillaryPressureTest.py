import OpenPNM
import scipy as sp


class CapillaryPressureTest:

    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        self.geo = OpenPNM.Geometry.GenericGeometry(network=self.net,
                                                    pores=self.net.Ps,
                                                    throats=self.net.Ts)
        self.geo['throat.diameter'] = 1
        self.geo['pore.diameter'] = 1
        self.water = OpenPNM.Phases.GenericPhase(network=self.net)
        self.phys = OpenPNM.Physics.GenericPhysics(network=self.net,
                                                   geometry=self.geo,
                                                   phase=self.water)

    def test_washburn_pore_values(self):
        self.water['pore.surface_tension'] = 0.072
        self.water['pore.contact_angle'] = 120
        f = OpenPNM.Physics.models.capillary_pressure.washburn
        self.phys.models.add(propname='pore.capillary_pressure',
                             model=f,
                             surface_tension='pore.surface_tension',
                             contact_angle='pore.contact_angle',
                             diameter='pore.diameter')
        self.phys.models.add(propname='throat.capillary_pressure',
                             model=f,
                             surface_tension='pore.surface_tension',
                             contact_angle='pore.contact_angle',
                             diameter='throat.diameter')
        a = 0.14399999999999993
        assert sp.allclose(self.water['pore.capillary_pressure'][0], a)
        self.phys.models.remove('pore.capillary_pressure')

    def test_washburn_throat_values(self):
        self.water['throat.surface_tension'] = 0.072
        self.water['throat.contact_angle'] = 120
        f = OpenPNM.Physics.models.capillary_pressure.washburn
        self.phys.models.add(propname='throat.capillary_pressure',
                             model=f,
                             surface_tension='throat.surface_tension',
                             contact_angle='throat.contact_angle',
                             diameter='throat.diameter')
        a = 0.14399999999999993
        assert sp.allclose(self.water['throat.capillary_pressure'][0], a)
        self.phys.models.remove('throat.capillary_pressure')

    def test_purcell_pore_values(self):
        self.water['pore.surface_tension'] = 0.072
        self.water['pore.contact_angle'] = 120
        f = OpenPNM.Physics.models.capillary_pressure.purcell
        self.phys.models.add(propname='pore.capillary_pressure',
                             model=f,
                             r_toroid=0.1,
                             surface_tension='pore.surface_tension',
                             contact_angle='pore.contact_angle',
                             diameter='pore.diameter')
        self.phys.models.add(propname='throat.capillary_pressure',
                             model=f,
                             r_toroid=0.1,
                             surface_tension='pore.surface_tension',
                             contact_angle='pore.contact_angle',
                             diameter='throat.diameter')
        a = 0.26206427646507374
        assert sp.allclose(self.water['pore.capillary_pressure'][0], a)
        self.phys.models.remove('pore.capillary_pressure')

    def test_purcell_throat_values(self):
        self.water['throat.surface_tension'] = 0.072
        self.water['throat.contact_angle'] = 120
        f = OpenPNM.Physics.models.capillary_pressure.purcell
        self.phys.models.add(propname='throat.capillary_pressure',
                             model=f,
                             r_toroid=0.1,
                             surface_tension='throat.surface_tension',
                             contact_angle='throat.contact_angle',
                             diameter='throat.diameter')
        a = 0.26206427646507374
        assert sp.allclose(self.water['throat.capillary_pressure'][0], a)
        self.phys.models.remove('throat.capillary_pressure')

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

    def test_cuboid(self):
        self.water['pore.surface_tension'] = 0.072
        self.water['pore.contact_angle'] = 120
        f = OpenPNM.Physics.models.capillary_pressure.cuboid
        self.phys.models.add(propname='throat.cuboid_pressure',
                             model=f)
        assert sp.sum(sp.isnan(self.phys['throat.cuboid_pressure'])) == 0

    def test_from_throat(self):
        f1 = OpenPNM.Physics.models.capillary_pressure.from_throat
        f2 = OpenPNM.Physics.models.capillary_pressure.washburn
        self.water['throat.surface_tension'] = 0.072
        self.water['throat.contact_angle'] = 120
        self.phys.models.add(propname='throat.capillary_pressure',
                             model=f2,
                             surface_tension='throat.surface_tension',
                             contact_angle='throat.contact_angle',
                             diameter='throat.diameter')
        self.phys.models.add(propname='pore.pc_min',
                             model=f1, operator='min')
        self.phys.models.add(propname='pore.pc_max',
                             model=f1, operator='max')
        self.phys.models.add(propname='pore.pc_mean',
                             model=f1, operator='mean')
        assert sp.all(self.phys['pore.pc_min'] <= self.phys['pore.pc_max'])
        assert sp.all((self.phys['pore.pc_mean'] <=
                      (self.phys['pore.pc_max']) *
                      (self.phys['pore.pc_mean'] >=
                       self.phys['pore.pc_min'])))

    def test_kelvin(self):
        f = OpenPNM.Physics.models.capillary_pressure.kelvin
        self.water['pore.temperature'] = 1.0
        self.water['pore.vapor_pressure'] = 1.0
        self.water['pore.molecular_weight'] = 1.0
        self.water['pore.density'] = 1.0
        self.water['pore.surface_tension'] = 1.0
        self.phys.models.add(propname='pore.pc_kelvin', model=f)
        assert sp.sum(sp.isnan(self.phys['pore.pc_kelvin'])) == 0
