import numpy as np
from numpy.testing import assert_approx_equal

import openpnm as op


class CapillaryPressureTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.net["throat.diameter"] = 1
        self.net["pore.diameter"] = 1
        self.water = op.phase.Phase(network=self.net)
        self.water["pore.surface_tension"] = 0.072
        self.water["pore.contact_angle"] = 120

    def test_washburn_pore_values(self):
        f = op.models.physics.capillary_pressure.washburn
        self.water.add_model(
            propname="pore.capillary_pressure",
            model=f,
            surface_tension="pore.surface_tension",
            contact_angle="pore.contact_angle",
            diameter="pore.diameter",
        )
        self.water.regenerate_models()
        a = 0.14399999999999993
        assert_approx_equal(self.water["pore.capillary_pressure"][0], a)
        del self.water["pore.capillary_pressure"]
        del self.water.models["pore.capillary_pressure@all"]

    def test_washburn_throat_values(self):
        f = op.models.physics.capillary_pressure.washburn
        self.water.add_model(
            propname="throat.capillary_pressure",
            model=f,
            surface_tension="throat.surface_tension",
            contact_angle="throat.contact_angle",
            diameter="throat.diameter",
        )
        self.water.regenerate_models()
        a = 0.14399999999999993
        assert_approx_equal(self.water["throat.capillary_pressure"][0], a)
        del self.water.models["throat.capillary_pressure"]

    def test_purcell_pore_values(self):
        f = op.models.physics.capillary_pressure.purcell
        self.water.add_model(
            propname="pore.capillary_pressure",
            model=f,
            r_toroid=0.1,
            surface_tension="pore.surface_tension",
            contact_angle="pore.contact_angle",
            diameter="pore.diameter",
        )
        self.water.regenerate_models()
        a = 0.2648436086476371
        assert_approx_equal(self.water["pore.capillary_pressure"][0], a)
        del self.water.models["pore.capillary_pressure"]

    def test_purcell_throat_values(self):
        f = op.models.physics.capillary_pressure.purcell
        self.water.add_model(
            propname="throat.capillary_pressure",
            model=f,
            r_toroid=0.1,
            surface_tension="throat.surface_tension",
            contact_angle="throat.contact_angle",
            diameter="throat.diameter",
        )
        self.water.regenerate_models()
        a = 0.2648436086476371
        assert_approx_equal(self.water["throat.capillary_pressure"][0], a)
        del self.water.models["throat.capillary_pressure"]

    # def test_purcell_bidirectional(self):
    #     f = op.models.physics.capillary_pressure.purcell_bidirectional
    #     self.net["pore.touch"] = (np.random.random(self.net.Np) + 0.5) * 0.1
    #     self.water.add_model(
    #         propname="throat.bidirectional",
    #         model=f,
    #         r_toroid=0.1,
    #         surface_tension="pore.surface_tension",
    #         contact_angle="pore.contact_angle",
    #         pore_diameter="pore.touch",
    #     )
    #     diff = (
    #         self.water["throat.bidirectional"][:, 0]
    #         - self.water["throat.bidirectional"][:, 1]
    #     )
    #     assert np.any(diff != 0)

    # def test_sinusoidal_bidirectional(self):
    #     f = op.models.physics.capillary_pressure.sinusoidal_bidirectional
    #     self.net["pore.touch"] = (np.random.random(self.net.Np) + 0.5) * 0.1
    #     self.net["throat.length"] = 1.0
    #     self.water.add_model(
    #         propname="throat.bidirectional",
    #         model=f,
    #         r_toroid=0.25,
    #         surface_tension="pore.surface_tension",
    #         contact_angle="pore.contact_angle",
    #         pore_diameter="pore.touch",
    #     )
    #     diff = (
    #         self.water["throat.bidirectional"][:, 0]
    #         - self.water["throat.bidirectional"][:, 1]
    #     )
    #     assert np.any(diff != 0)

    # def test_ransohoff_snapoff_verts(self):
    #     ws = op.Workspace()
    #     ws.clear()
    #     bp = np.array(
    #         [
    #             [0.25, 0.25, 0.25],
    #             [0.25, 0.75, 0.25],
    #             [0.75, 0.25, 0.25],
    #             [0.75, 0.75, 0.25],
    #             [0.25, 0.25, 0.75],
    #             [0.25, 0.75, 0.75],
    #             [0.75, 0.25, 0.75],
    #             [0.75, 0.75, 0.75],
    #         ]
    #     )
    #     scale = 1e-4
    #     np.random.seed(1)
    #     p = (np.random.random([len(bp), 3]) - 0.5) / 1000
    #     bp += p
    #     fiber_rad = 2e-6
    #     bp = op.topotools.reflect_base_points(points=bp, domain_size=[1, 1, 1])
    #     prj = op.materials.VoronoiFibers(
    #         fiber_rad=fiber_rad,
    #         resolution=1e-6,
    #         shape=[scale, scale, scale],
    #         points=bp * scale,
    #         name="test",
    #     )
    #     net = prj.network
    #     del_geom = prj.geometries()["test_del"]
    #     vor_geom = prj.geometries()["test_vor"]
    #     f = op.models.physics.capillary_pressure.ransohoff_snap_off
    #     water = op.phase.Phase(network=net)
    #     water["pore.surface_tension"] = 0.072
    #     water["pore.contact_angle"] = 45
    #     phys1 = op.physics.GenericPhysics(network=net, geometry=del_geom, phase=water)
    #     phys1.add_model(propname="throat.snap_off", model=f, wavelength=fiber_rad)
    #     phys1.add_model(
    #         propname="throat.snap_off_pair",
    #         model=f,
    #         wavelength=fiber_rad,
    #         require_pair=True,
    #     )
    #     phys2 = op.physics.GenericPhysics(network=net, geometry=vor_geom, phase=water)
    #     phys2.add_model(propname="throat.snap_off", model=f, wavelength=fiber_rad)
    #     phys2.add_model(
    #         propname="throat.snap_off_pair",
    #         model=f,
    #         wavelength=fiber_rad,
    #         require_pair=True,
    #     )
    #     ts = ~net["throat.interconnect"]
    #     assert ~np.any(np.isnan(water["throat.snap_off"][ts]))
    #     assert np.any(~np.isnan(water["throat.snap_off_pair"][ts]))


if __name__ == "__main__":

    t = CapillaryPressureTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith("test"):
            print(f"Running test: {item}")
            t.__getattribute__(item)()
