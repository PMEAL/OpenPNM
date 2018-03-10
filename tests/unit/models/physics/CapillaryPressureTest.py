import openpnm as op
import scipy as sp


class CapillaryPressureTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.geo = op.geometry.GenericGeometry(network=self.net,
                                               pores=self.net.Ps,
                                               throats=self.net.Ts)
        self.geo['throat.diameter'] = 1
        self.geo['pore.diameter'] = 1
        self.water = op.phases.GenericPhase(network=self.net)
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              geometry=self.geo,
                                              phase=self.water)

    def test_washburn_pore_values(self):
        self.water['pore.surface_tension'] = 0.072
        self.water['pore.contact_angle'] = 120
        f = op.models.physics.capillary_pressure.washburn
        self.phys.add_model(propname='pore.capillary_pressure',
                            model=f,
                            surface_tension='pore.surface_tension',
                            contact_angle='pore.contact_angle',
                            diameter='pore.diameter')
        self.phys.add_model(propname='throat.capillary_pressure',
                            model=f,
                            surface_tension='pore.surface_tension',
                            contact_angle='pore.contact_angle',
                            diameter='throat.diameter')
        self.phys.regenerate_models()
        a = 0.14399999999999993
        assert sp.allclose(self.water['pore.capillary_pressure'][0], a)
        self.phys.models.remove('pore.capillary_pressure')

    def test_washburn_throat_values(self):
        self.water['throat.surface_tension'] = 0.072
        self.water['throat.contact_angle'] = 120
        f = op.models.physics.capillary_pressure.washburn
        self.phys.add_model(propname='throat.capillary_pressure',
                            model=f,
                            surface_tension='throat.surface_tension',
                            contact_angle='throat.contact_angle',
                            diameter='throat.diameter')
        self.phys.regenerate_models()
        a = 0.14399999999999993
        assert sp.allclose(self.water['throat.capillary_pressure'][0], a)
        self.phys.models.remove('throat.capillary_pressure')

    def test_purcell_pore_values(self):
        self.water['pore.surface_tension'] = 0.072
        self.water['pore.contact_angle'] = 120
        f = op.models.physics.capillary_pressure.purcell
        self.phys.add_model(propname='pore.capillary_pressure',
                            model=f,
                            r_toroid=0.1,
                            surface_tension='pore.surface_tension',
                            contact_angle='pore.contact_angle',
                            diameter='pore.diameter')
        self.phys.add_model(propname='throat.capillary_pressure',
                            model=f,
                            r_toroid=0.1,
                            surface_tension='pore.surface_tension',
                            contact_angle='pore.contact_angle',
                            diameter='throat.diameter')
        self.phys.regenerate_models()
        a = 0.26206427646507374
        assert sp.allclose(self.water['pore.capillary_pressure'][0], a)
        self.phys.models.remove('pore.capillary_pressure')

    def test_purcell_throat_values(self):
        self.water['throat.surface_tension'] = 0.072
        self.water['throat.contact_angle'] = 120
        f = op.models.physics.capillary_pressure.purcell
        self.phys.add_model(propname='throat.capillary_pressure',
                            model=f,
                            r_toroid=0.1,
                            surface_tension='throat.surface_tension',
                            contact_angle='throat.contact_angle',
                            diameter='throat.diameter')
        self.phys.regenerate_models()
        a = 0.26206427646507374
        assert sp.allclose(self.water['throat.capillary_pressure'][0], a)
        self.phys.models.remove('throat.capillary_pressure')

    def test_static_pressure(self):
        self.water['pore.density'] = 1000
        self.water['pore.occupancy'] = True
        f = op.models.physics.capillary_pressure.static_pressure
        self.phys.add_model(propname='pore.static_pressure',
                            model=f,
                            gravity=[0, 0, 9.81],
                            pore_density='pore.density',
                            pore_occupancy='pore.occupancy')
        self.phys.regenerate_models()
        a = [0., 9810., 19620.]
        assert sp.all(sp.unique(self.water['pore.static_pressure']) == a)

    def test_cuboid(self):
        self.water['pore.surface_tension'] = 0.072
        self.water['pore.contact_angle'] = 120
        f = op.models.physics.capillary_pressure.cuboid
        self.phys.add_model(propname='throat.cuboid_pressure',
                            model=f)
        self.phys.regenerate_models()
        assert sp.sum(sp.isnan(self.phys['throat.cuboid_pressure'])) == 0

    def test_from_throat(self):
        f1 = op.models.physics.capillary_pressure.from_throat
        f2 = op.models.physics.capillary_pressure.washburn
        self.water['throat.surface_tension'] = 0.072
        self.water['throat.contact_angle'] = 120
        self.phys.add_model(propname='throat.capillary_pressure',
                            model=f2,
                            surface_tension='throat.surface_tension',
                            contact_angle='throat.contact_angle',
                            diameter='throat.diameter')
        self.phys.add_model(propname='pore.pc_min',
                            model=f1, operator='min')
        self.phys.add_model(propname='pore.pc_max',
                            model=f1, operator='max')
        self.phys.add_model(propname='pore.pc_mean',
                            model=f1, operator='mean')
        self.phys.regenerate_models()
        assert sp.all(self.phys['pore.pc_min'] <= self.phys['pore.pc_max'])
        assert sp.all((self.phys['pore.pc_mean'] <=
                      (self.phys['pore.pc_max']) *
                      (self.phys['pore.pc_mean'] >=
                       self.phys['pore.pc_min'])))

    def test_kelvin(self):
        f = op.models.physics.capillary_pressure.kelvin
        self.water['pore.temperature'] = 1.0
        self.water['pore.vapor_pressure'] = 1.0
        self.water['pore.molecular_weight'] = 1.0
        self.water['pore.density'] = 1.0
        self.water['pore.surface_tension'] = 1.0
        self.phys.add_model(propname='pore.pc_kelvin', model=f)
        assert sp.sum(sp.isnan(self.phys['pore.pc_kelvin'])) == 0


if __name__ == '__main__':

    t = CapillaryPressureTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
