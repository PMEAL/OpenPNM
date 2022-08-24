import pytest
import numpy as np
from numpy.testing import assert_allclose
import openpnm as op
import openpnm.models.geometry.conduit_lengths as mods
import openpnm.models.geometry as gm


class ConduitLengthsTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 1, 1], spacing=1.0)
        self.net.regenerate_models()
        self.air = op.phase.Air(network=self.net)
        self.net["throat.diameter"] = 0.4
        self.net["pore.diameter"] = [0.5, 0.9, 0.7]

    def test_spheres_and_cylinders(self):
        L_actual = mods.spheres_and_cylinders(self.net)
        L_desired = np.array([[0.15      , 0.44688711, 0.40311289],
                              [0.40311289, 0.30965898, 0.28722813]])
        assert_allclose(L_actual, L_desired)
        # Incompatible data with model assumptions
        self.net["pore.diameter"][1] = 0.3
        with pytest.raises(Exception):
            L_actual = mods.spheres_and_cylinders(self.net)
        self.net["pore.diameter"][1] = 0.9
        # Overlapping pores
        self.net["pore.diameter"][0] = 1.2
        L_actual = mods.spheres_and_cylinders(self.net)
        L_desired = np.array([[5.78750000e-01, 1.00000000e-15, 4.21250000e-01],
                              [4.03112887e-01, 3.09658980e-01, 2.87228132e-01]])
        assert_allclose(L_actual, L_desired)
        self.net["pore.diameter"][0] = 0.5

    def test_circles_and_rectangles(self):
        L_actual = mods.circles_and_rectangles(self.net)
        L_desired = np.array([[0.15      , 0.44688711, 0.40311289],
                              [0.40311289, 0.30965898, 0.28722813]])
        assert_allclose(L_actual, L_desired)
        # Incompatible data with model assumptions
        self.net["pore.diameter"][1] = 0.3
        with pytest.raises(Exception):
            L_actual = mods.circles_and_squares(self.net)
        self.net["pore.diameter"][1] = 0.9
        # Overlapping pores
        self.net["pore.diameter"][0] = 1.2
        L_actual = mods.circles_and_rectangles(self.net)
        L_desired = np.array([[5.78750000e-01, 1.00000000e-15, 4.21250000e-01],
                              [4.03112887e-01, 3.09658980e-01, 2.87228132e-01]])
        assert_allclose(L_actual, L_desired)
        self.net["pore.diameter"][0] = 0.5

    def test_cones_and_cylinders(self):
        L_actual = mods.cones_and_cylinders(self.net)
        L_desired = np.array([[0.25, 0.3, 0.45],
                              [0.45, 0.2, 0.35]])
        assert_allclose(L_actual, L_desired)
        # Overlapping pores
        self.net["pore.diameter"][0] = 1.2
        L_actual = mods.cones_and_cylinders(self.net)
        L_desired = np.array([[5.7875e-01, 1.0000e-15, 4.2125e-01],
                              [4.5000e-01, 2.0000e-01, 3.5000e-01]])
        assert_allclose(L_actual, L_desired)
        self.net["pore.diameter"][0] = 0.5

    def test_trapezoids_and_rectangles(self):
        L_actual = mods.trapezoids_and_rectangles(self.net)
        L_desired = np.array([[0.25, 0.3, 0.45],
                              [0.45, 0.2, 0.35]])
        assert_allclose(L_actual, L_desired)
        # Overlapping pores
        self.net["pore.diameter"][0] = 1.2
        L_actual = mods.trapezoids_and_rectangles(self.net)
        L_desired = np.array([[5.7875e-01, 1.0000e-15, 4.2125e-01],
                              [4.5000e-01, 2.0000e-01, 3.5000e-01]])
        assert_allclose(L_actual, L_desired)
        self.net["pore.diameter"][0] = 0.5

    def test_pyramids_and_cuboids(self):
        L_actual = mods.pyramids_and_cuboids(self.net)
        L_desired = np.array([[0.25, 0.3, 0.45],
                              [0.45, 0.2, 0.35]])
        assert_allclose(L_actual, L_desired)
        # Overlapping pores
        self.net["pore.diameter"][0] = 1.2
        L_actual = mods.pyramids_and_cuboids(self.net)
        L_desired = np.array([[5.7875e-01, 1.0000e-15, 4.2125e-01],
                              [4.5000e-01, 2.0000e-01, 3.5000e-01]])
        assert_allclose(L_actual, L_desired)
        self.net["pore.diameter"][0] = 0.5

    def test_cubes_and_cuboids(self):
        L_actual = mods.cubes_and_cuboids(self.net)
        L_desired = np.array([[0.25, 0.3, 0.45],
                              [0.45, 0.2, 0.35]])
        assert_allclose(L_actual, L_desired)
        # Overlapping pores
        self.net["pore.diameter"][0] = 1.2
        L_actual = mods.cubes_and_cuboids(self.net)
        L_desired = np.array([[6.0e-01, 1.0e-15, 4.0e-01],
                              [4.5e-01, 2.0e-01, 3.5e-01]])
        assert_allclose(L_actual, L_desired)
        self.net["pore.diameter"][0] = 0.5

    def test_squares_and_rectangles(self):
        L_actual = mods.squares_and_rectangles(self.net)
        L_desired = np.array([[0.25, 0.3, 0.45],
                              [0.45, 0.2, 0.35]])
        assert_allclose(L_actual, L_desired)
        # Overlapping pores
        self.net["pore.diameter"][0] = 1.2
        L_actual = mods.cubes_and_cuboids(self.net)
        L_desired = np.array([[6.0e-01, 1.0e-15, 4.0e-01],
                              [4.5e-01, 2.0e-01, 3.5e-01]])
        assert_allclose(L_actual, L_desired)
        self.net["pore.diameter"][0] = 0.5

    def test_intersecting_cones(self):
        self.net.add_model(propname='throat.coords',
                           model=gm.throat_centroid.pore_coords)
        L_actual = mods.intersecting_cones(self.net)
        L_desired = np.array([[0.5, 0., 0.5],
                              [0.5, 0., 0.5]])
        assert_allclose(L_actual, L_desired)

    def test_hybrid_cones_and_cylinders(self):
        self.net.add_model(propname='throat.coords',
                           model=gm.throat_centroid.pore_coords)
        self.net["pore.diameter"][0] = 1.2
        L_actual = mods.hybrid_cones_and_cylinders(self.net)
        L_desired = np.array([[0.5, 0., 0.5],
                              [0.45, 0.2, 0.35]])
        assert_allclose(L_actual, L_desired)
        self.net["pore.diameter"][0] = 0.5

    def test_intersecting_pyramids(self):
        self.net.add_model(propname='throat.coords',
                           model=gm.throat_centroid.pore_coords)
        L_actual = mods.intersecting_pyramids(self.net)
        L_desired = np.array([[0.5, 0., 0.5],
                              [0.5, 0., 0.5]])
        assert_allclose(L_actual, L_desired)

    def test_hybrid_pyramids_and_cuboids(self):
        self.net.add_model(propname='throat.coords',
                           model=gm.throat_centroid.pore_coords)
        self.net["pore.diameter"][0] = 1.2
        L_actual = mods.hybrid_pyramids_and_cuboids(self.net)
        L_desired = np.array([[0.5, 0., 0.5],
                              [0.45, 0.2, 0.35]])
        assert_allclose(L_actual, L_desired)
        self.net["pore.diameter"][0] = 0.5

    def test_intersecting_trapezoids(self):
        self.net.add_model(propname='throat.coords',
                           model=gm.throat_centroid.pore_coords)
        L_actual = mods.intersecting_trapezoids(self.net)
        L_desired = np.array([[0.5, 0., 0.5],
                              [0.5, 0., 0.5]])
        assert_allclose(L_actual, L_desired)

    def test_hybrid_trapezoids_and_rectangles(self):
        self.net.add_model(propname='throat.coords',
                           model=gm.throat_centroid.pore_coords)
        self.net["pore.diameter"][0] = 1.2
        L_actual = mods.hybrid_trapezoids_and_rectangles(self.net)
        L_desired = np.array([[0.5, 0., 0.5],
                              [0.45, 0.2, 0.35]])
        assert_allclose(L_actual, L_desired)
        self.net["pore.diameter"][0] = 0.5


if __name__ == '__main__':

    t = ConduitLengthsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
