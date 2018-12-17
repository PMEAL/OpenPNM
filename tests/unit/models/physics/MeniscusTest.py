import openpnm as op
import openpnm.models.physics as pm
import scipy as sp


class MeniscusTest:

    def setup_class(self):
        sp.random.seed(1)
        self.net = op.network.Cubic(shape=[5, 1, 5], spacing=5e-5)
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=self.net.pores(),
                                            throats=self.net.throats())
        self.phase = op.phases.Water(network=self.net)
        self.phys = op.physics.Standard(network=self.net,
                                        phase=self.phase,
                                        geometry=self.geo)

    def test_toroidal_touch(self):
        phys = self.phys
        r_tor = 1e-6
        self.geo['throat.touch_length'] = 2e-6
        phys.add_model(propname='throat.tor_max',
                       model=pm.meniscus.purcell,
                       mode='max',
                       r_toroid=r_tor)
        phys.add_model(propname='throat.tor_touch',
                       model=pm.meniscus.purcell,
                       mode='touch',
                       r_toroid=r_tor)
        assert sp.any(phys['throat.tor_touch'] < phys['throat.tor_max'])

    def test_sinusoidal_touch(self):
        phys = self.phys
        self.geo['throat.amplitude'] = 5e-6
        self.geo['throat.touch_length'] = 1e-6
        phys.add_model(propname='throat.sin_pressure_max',
                       model=pm.meniscus.sinusoidal,
                       mode='max')
        phys.add_model(propname='throat.sin_pressure_touch',
                       model=pm.meniscus.sinusoidal,
                       mode='touch')
        h = phys.check_data_health()
        for check in h.values():
            if len(check) > 0:
                assert 1 == 2
        assert sp.any((phys['throat.sin_pressure_touch'] <
                       phys['throat.sin_pressure_max']))

    def test_sinusoidal(self):
        phys = self.phys
        self.geo['throat.amplitude'] = 5e-6
        phys.add_model(propname='throat.sin_pressure',
                       model=pm.meniscus.sinusoidal,
                       mode='max')
        phys.add_model(propname='throat.sin_meniscus',
                       model=pm.meniscus.sinusoidal,
                       mode='men',
                       target_Pc=5000)
        h = phys.check_data_health()
        for check in h.values():
            if len(check) > 0:
                assert 1 == 2

    def test_toroidal(self):
        phys = self.phys
        r_tor = 1e-6
        phys.add_model(propname='throat.purcell_pressure',
                       model=pm.capillary_pressure.purcell,
                       r_toroid=r_tor)
        phys.add_model(propname='throat.tor_pressure',
                       model=pm.meniscus.purcell,
                       mode='max',
                       r_toroid=r_tor,
                       num_points=1000)
        phys.add_model(propname='throat.tor_meniscus',
                       model=pm.meniscus.purcell,
                       mode='men',
                       r_toroid=r_tor,
                       target_Pc=5000)
        a = sp.around(phys['throat.purcell_pressure'], 10)
        b = sp.around(phys['throat.tor_pressure'], 10)
        assert sp.allclose(a, b)
        h = phys.check_data_health()
        for check in h.values():
            if len(check) > 0:
                assert 1 == 2

    def test_general_toroidal(self):
        phys = self.phys
        r_tor = 1e-6
        phys.add_model(propname='throat.purcell_pressure',
                       model=pm.capillary_pressure.purcell,
                       r_toroid=r_tor)
        phys['throat.scale_a'] = r_tor
        phys['throat.scale_b'] = r_tor
        phys.add_model(propname='throat.general_pressure',
                       model=pm.meniscus.general_toroidal,
                       mode='max',
                       num_points=1000)
        a = sp.around(phys['throat.purcell_pressure'], 10)
        b = sp.around(phys['throat.general_pressure'], 10)
        assert sp.allclose(a, b)
        h = phys.check_data_health()
        for check in h.values():
            if len(check) > 0:
                assert 1 == 2


if __name__ == '__main__':

    t = MeniscusTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
