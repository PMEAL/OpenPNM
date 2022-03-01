import numpy as np
import openpnm as op
import openpnm.models.physics as pm


class MeniscusTest:
    def setup_class(self):
        np.random.seed(1)
        self.net = op.network.Cubic(shape=[5, 1, 5], spacing=5e-5)
        self.geo = op.geometry.SpheresAndCylinders(network=self.net,
                                                   pores=self.net.pores(),
                                                   throats=self.net.throats())
        self.phase = op.phase.Water(network=self.net)
        self.phys = op.physics.Standard(network=self.net,
                                        phase=self.phase,
                                        geometry=self.geo)

    def test_toroidal_touch(self):
        phys = self.phys
        r_tor = 1e-6
        self.geo["throat.touch_length"] = 2e-6
        phys.add_model(propname="throat.tor_max",
                       model=pm.meniscus.purcell,
                       mode="max",
                       r_toroid=r_tor,)
        phys.add_model(propname="throat.tor_touch",
                           model=pm.meniscus.purcell,
                           mode="touch",
                           r_toroid=r_tor,)
        assert np.any(phys["throat.tor_touch"] < phys["throat.tor_max"])

    def test_sinusoidal_touch(self):
        phys = self.phys
        self.geo["throat.amplitude"] = 5e-6
        self.geo["throat.touch_length"] = 1e-6
        phys.add_model(propname="throat.sin_pressure_max",
                       model=pm.meniscus.sinusoidal, mode="max")
        phys.add_model(propname="throat.sin_pressure_touch",
                       model=pm.meniscus.sinusoidal,
                       mode="touch")
        h = phys.project.check_data_health(phys)
        for check in h.values():
            if len(check) > 0:
                assert 1 == 2
        assert np.any(
            (phys["throat.sin_pressure_touch"] < phys["throat.sin_pressure_max"])
        )

    def test_sinusoidal(self):
        phys = self.phys
        self.geo["throat.amplitude"] = 5e-6
        phys.add_model(propname="throat.sin_pressure",
                       model=pm.meniscus.sinusoidal,
                       mode="max")
        phys.add_model(propname="throat.sin_meniscus",
                       model=pm.meniscus.sinusoidal,
                       mode="men",
                       target_Pc=5000)
        h = phys.project.check_data_health(phys)
        for check in h.values():
            if len(check) > 0:
                assert 1 == 2

    def test_toroidal(self):
        phys = self.phys
        r_tor = 1e-6
        phys.add_model(propname="throat.purcell_pressure",
                       model=pm.capillary_pressure.purcell,
                       r_toroid=r_tor)
        phys.add_model(propname="throat.tor_pressure",
                       model=pm.meniscus.purcell,
                       mode="max",
                       r_toroid=r_tor,
                       num_points=1000)
        phys.add_model(propname="throat.tor_meniscus",
                       model=pm.meniscus.purcell,
                       mode="men",
                       r_toroid=r_tor,
                       target_Pc=5000)
        a = np.around(phys["throat.purcell_pressure"], 10)
        b = np.around(phys["throat.tor_pressure"], 10)
        assert np.allclose(a, b)
        h = phys.project.check_data_health(phys)
        for check in h.values():
            if len(check) > 0:
                assert 1 == 2

    def test_general_toroidal(self):
        phys = self.phys
        r_tor = 1e-6
        phys.add_model(propname="throat.purcell_pressure",
                       model=pm.capillary_pressure.purcell,
                       r_toroid=r_tor)
        phys["throat.scale_a"] = r_tor
        phys["throat.scale_b"] = r_tor
        phys.add_model(propname="throat.general_pressure",
                       model=pm.meniscus.general_toroidal,
                       mode="max",
                       num_points=1000)
        a = np.around(phys["throat.purcell_pressure"], 10)
        b = np.around(phys["throat.general_pressure"], 10)
        assert np.allclose(a, b)
        h = phys.project.check_data_health(phys)
        for check in h.values():
            if len(check) > 0:
                assert 1 == 2

    def test_exceptions(self):
        phys = self.phys
        r_tor = 1e-6
        phys["throat.scale_a"] = r_tor
        phys["throat.scale_b"] = r_tor
        phys.add_model(propname="throat.elliptical_pressure",
                       model=pm.meniscus.general_toroidal,
                       mode="max",
                       profile_equation="elliptical",
                       num_points=1000)
        phys.add_model(propname="throat.exception_pressure",
                       model=pm.meniscus.general_toroidal,
                       mode="max",
                       profile_equation="scooby-doo",
                       num_points=1000)
        a = np.around(phys["throat.elliptical_pressure"], 10)
        b = np.around(phys["throat.exception_pressure"], 10)
        assert np.allclose(a, b)
        phys.add_model(propname="throat.no_target_pressure",
                       model=pm.meniscus.general_toroidal,
                       mode="men",
                       num_points=1000)
        phys.add_model(propname="throat.small_target_pressure",
                       model=pm.meniscus.general_toroidal,
                       mode="men",
                       target_Pc=1.0e-7,
                       num_points=1000)
        a = np.around(phys["throat.no_target_pressure.radius"], 10)
        b = np.around(phys["throat.small_target_pressure.radius"], 10)
        assert np.allclose(a, b)
        h = phys.project.check_data_health(phys)
        for check in h.values():
            if len(check) > 0:
                assert 1 == 2


if __name__ == "__main__":

    t = MeniscusTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith("test"):
            print("running test: " + item)
            t.__getattribute__(item)()
