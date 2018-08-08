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

    def test_sinusoidal(self):
        phys = self.phys
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

    def test_purcell(self):
        phys = self.phys
        r_tor = 1e-6
        phys.add_model(propname='throat.purcell_pressure',
                       model=pm.capillary_pressure.purcell,
                       r_toroid=r_tor)
        phys.add_model(propname='throat.purcell_men_alpha',
                       model=pm.meniscus.purcell_filling_angle,
                       r_toroid=r_tor,
                       Pc=5000)
        phys.add_model(propname='throat.purcell_men_rad',
                       model=pm.meniscus.purcell_meniscus_radius,
                       filling_angle='throat.purcell_men_alpha',
                       r_toroid=r_tor)
        phys.add_model(propname='throat.purcell_men_cen',
                       model=pm.meniscus.purcell_meniscus_center,
                       filling_angle='throat.purcell_men_alpha',
                       men_rad='throat.purcell_men_rad',
                       r_toroid=r_tor)
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
