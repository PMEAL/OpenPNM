from numpy import array
from numpy.testing import assert_allclose
import openpnm as op
import openpnm.models.geometry.hydraulic_size_factors as mods


class ConduitHydraulicSizeFactorsTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 1, 1], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(
            network=self.net, pores=self.net.Ps, throats=self.net.Ts
        )
        self.air = op.phase.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(
            network=self.net, phase=self.air, geometry=self.geo
        )
        self.geo["throat.diameter"] = 0.4
        self.geo["pore.diameter"] = [1.2, 0.9, 0.7]

    def test_spheres_and_cylinders(self):
        S_actual = mods.spheres_and_cylinders(self.geo)
        S_desired = {
            'pore1': array([0.01068901, 0.01195694]),
            'throat': array([6.28318531e+11, 2.02906607e-03]),
            'pore2': array([0.00771764, 0.00917032])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k], rtol=1e-5)

    def test_circles_and_rectangles(self):
        S_actual = mods.circles_and_rectangles(self.geo)
        S_desired = {
            'pore1': array([0.06563123, 0.06697876]),
            'throat': array([5.33333333e+12, 1.72232477e-02]),
            'pore2': array([0.05072058, 0.05686537])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k], rtol=1e-5)

    def test_cones_and_cylinders(self):
        S_actual = mods.cones_and_cylinders(self.geo)
        S_desired = {
            'pore1': array([0.00676442, 0.0057399]),
            'throat': array([6.28318531e+11, 3.14159265e-03]),
            'pore2': array([0.00613165, 0.00496574])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k], rtol=1e-5)

    def test_trapezoids_and_rectangles(self):
        S_actual = mods.trapezoids_and_rectangles(self.geo)
        S_desired = {
            'pore1': array([0.04146868, 0.03692308]),
            'throat': array([5.33333333e+12, 2.66666667e-02]),
            'pore2': array([0.03944305, 0.03393939])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k], rtol=1e-5)

    def test_pyramids_and_cuboids(self):
        S_actual = mods.pyramids_and_cuboids(self.geo)
        S_desired = {
            'pore1': array([0.01047182, 0.00888579]),
            'throat': array([9.72683363e+11, 4.86341681e-03]),
            'pore2': array([0.00949224, 0.00768734])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k], rtol=1e-5)

    def test_cubes_and_cuboids(self):
        S_actual = mods.cubes_and_cuboids(self.geo)
        S_desired = {
            'pore1': array([0.13131225, 0.05539736]),
            'throat': array([9.72683363e+11, 4.86341681e-03]),
            'pore2': array([0.06232203, 0.02606487])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k], rtol=1e-5)

    def test_squares_and_rectangles(self):
        S_actual = mods.squares_and_rectangles(self.geo)
        S_desired = {
            'pore1': array([0.24 , 0.135]),
            'throat': array([5.33333333e+12, 2.66666667e-02]),
            'pore2': array([0.151875  , 0.08166667])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k], rtol=1e-5)

    def test_ncylinders_in_series(self):
        S_actual = mods.ncylinders_in_series(self.geo)
        S_desired = {
            'pore1': array([0.0020471 , 0.00612218]),
            'throat': array([6.28318531e+11, 2.02906607e-03]),
            'pore2': array([0.00261812, 0.00661558])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k], rtol=1e-5)


if __name__ == '__main__':

    t = ConduitHydraulicSizeFactorsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
