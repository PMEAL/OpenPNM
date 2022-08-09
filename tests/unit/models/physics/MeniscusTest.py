import numpy as np
import openpnm as op
import openpnm.models.physics as pm


class MeniscusTest:
    def setup_class(self):
        np.random.seed(1)
        self.net = op.network.Cubic(shape=[5, 1, 5], spacing=5e-5)
        self.net.add_model_collection(
            op.models.collections.geometry.spheres_and_cylinders)
        self.net.regenerate_models()
        self.phase = op.phase.Water(network=self.net)

    def test_toroidal_touch(self):
        r_tor = 1e-6
        self.net["throat.touch_length"] = 2e-6
        self.phase.add_model(propname="throat.tor_max",
                             model=pm.meniscus.purcell,
                             mode="max",
                             r_toroid=r_tor,)
        self.phase.add_model(propname="throat.tor_touch",
                             model=pm.meniscus.purcell,
                             mode="touch",
                             r_toroid=r_tor,)
        assert np.any(self.phase["throat.tor_touch"] < self.phase["throat.tor_max"])

    def test_sinusoidal_touch(self):
        self.net["throat.amplitude"] = 5e-6
        self.net["throat.touch_length"] = 1e-6
        self.phase.add_model(propname="throat.sin_pressure_max",
                             model=pm.meniscus.sinusoidal, mode="max")
        self.phase.add_model(propname="throat.sin_pressure_touch",
                             model=pm.meniscus.sinusoidal,
                             mode="touch")
        assert np.any(
            (self.phase["throat.sin_pressure_touch"]
              < self.phase["throat.sin_pressure_max"])
        )
        # h = self.project.check_data_health(phys)
        # for check in h.values():
            # if len(check) > 0:
                # assert 1 == 2

    def test_sinusoidal(self):
        self.net["throat.amplitude"] = 5e-6
        self.phase.add_model(propname="throat.sin_pressure",
                             model=pm.meniscus.sinusoidal,
                             mode="max")
        self.phase.add_model(propname="throat.sin_meniscus",
                             model=pm.meniscus.sinusoidal,
                             mode="men",
                             target_Pc=5000)
        # h = self.phase.project.check_data_health(phys)
        # for check in h.values():
        #     if len(check) > 0:
        #         assert 1 == 2

    def test_toroidal(self):
        r_tor = 1e-6
        self.phase.add_model(propname="throat.purcell_pressure",
                             model=pm.capillary_pressure.purcell,
                             r_toroid=r_tor)
        self.phase.add_model(propname="throat.tor_pressure",
                             model=pm.meniscus.purcell,
                             mode="max",
                             r_toroid=r_tor,
                             num_points=1000)
        self.phase.add_model(propname="throat.tor_meniscus",
                             model=pm.meniscus.purcell,
                             mode="men",
                             r_toroid=r_tor,
                             target_Pc=5000)
        a = np.around(self.phase["throat.purcell_pressure"], 10)
        b = np.around(self.phase["throat.tor_pressure"], 10)
        assert np.allclose(a, b)
        # h = phys.project.check_data_health(phys)
        # for check in h.values():
        #     if len(check) > 0:
        #         assert 1 == 2

    def test_general_toroidal(self):
        r_tor = 1e-6
        self.phase.add_model(propname="throat.purcell_pressure",
                             model=pm.capillary_pressure.purcell,
                             r_toroid=r_tor)
        self.phase["throat.scale_a"] = r_tor
        self.phase["throat.scale_b"] = r_tor
        self.phase.add_model(propname="throat.general_pressure",
                             model=pm.meniscus.general_toroidal,
                             mode="max",
                             num_points=1000)
        a = np.around(self.phase["throat.purcell_pressure"], 10)
        b = np.around(self.phase["throat.general_pressure"], 10)
        assert np.allclose(a, b)
        # h = phys.project.check_data_health(phys)
        # for check in h.values():
        #     if len(check) > 0:
        #         assert 1 == 2

    def test_exceptions(self):
        r_tor = 1e-6
        self.phase["throat.scale_a"] = r_tor
        self.phase["throat.scale_b"] = r_tor
        self.phase.add_model(propname="throat.elliptical_pressure",
                             model=pm.meniscus.general_toroidal,
                             mode="max",
                             profile_equation="elliptical",
                             num_points=1000)
        self.phase.add_model(propname="throat.exception_pressure",
                             model=pm.meniscus.general_toroidal,
                             mode="max",
                             profile_equation="scooby-doo",
                             num_points=1000)
        a = np.around(self.phase["throat.elliptical_pressure"], 10)
        b = np.around(self.phase["throat.exception_pressure"], 10)
        assert np.allclose(a, b)
        self.phase.add_model(propname="throat.no_target_pressure",
                             model=pm.meniscus.general_toroidal,
                             mode="men",
                             num_points=1000)
        self.phase.add_model(propname="throat.small_target_pressure",
                             model=pm.meniscus.general_toroidal,
                             mode="men",
                             target_Pc=1.0e-7,
                             num_points=1000)
        a = np.around(self.phase["throat.no_target_pressure.radius"], 10)
        b = np.around(self.phase["throat.small_target_pressure.radius"], 10)
        # assert np.allclose(a, b)
        # h = phys.project.check_data_health(phys)
        # for check in h.values():
        #     if len(check) > 0:
        #         assert 1 == 2


if __name__ == "__main__":

    t = MeniscusTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith("test"):
            print("running test: " + item)
            t.__getattribute__(item)()
