import openpnm as op
from numpy.testing import assert_approx_equal


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
        assert_approx_equal(self.water['pore.capillary_pressure'][0], a)
        self.phys.remove_model('pore.capillary_pressure')

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
        assert_approx_equal(self.water['throat.capillary_pressure'][0], a)
        self.phys.remove_model('throat.capillary_pressure')

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
        assert_approx_equal(self.water['pore.capillary_pressure'][0], a)
        self.phys.remove_model('pore.capillary_pressure')

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
        assert_approx_equal(self.water['throat.capillary_pressure'][0], a)
        self.phys.remove_model('throat.capillary_pressure')


if __name__ == '__main__':

    t = CapillaryPressureTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
