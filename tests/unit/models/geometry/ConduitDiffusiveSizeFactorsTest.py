from numpy import array
from numpy.testing import assert_allclose
import openpnm as op
import openpnm.models.geometry.diffusive_size_factors as mods


class ConduitDiffusiveSizeFactorsTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 1, 1], spacing=1.0)
        self.geo = op.geometry.GenericGeometry(
            network=self.net, pores=self.net.Ps, throats=self.net.Ts
        )
        self.air = op.phases.Air(network=self.net)
        self.phys = op.physics.GenericPhysics(
            network=self.net, phase=self.air, geometry=self.geo
        )
        self.geo["throat.diameter"] = 0.4
        self.geo["pore.diameter"] = [1.2, 0.9, 0.7]

    def test_spheres_and_cylinders(self):
        S_actual = mods.spheres_and_cylinders(self.geo)
        S_desired = {
            'pore1': array([0.93875728, 0.97459088]),
            'throat': array([1.25663706e+14, 4.05813214e-01]),
            'pore2': array([0.82884551, 0.94886745])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k])

    def test_circles_and_rectangles(self):
        S_actual = mods.circles_and_rectangles(self.geo)
        S_desired = {
            'pore1': array([1.53390798, 1.80140852]),
            'throat': array([4.00000000e+14, 1.29174358e+00]),
            'pore2': array([1.65097536, 2.07781253])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k])

    def test_cones_and_cylinders(self):
        S_actual = mods.cones_and_cylinders(self.geo)
        S_desired = {
            'pore1': array([0.65138854, 0.62831853]),
            'throat': array([1.25663706e+14, 6.28318531e-01]),
            'pore2': array([0.6712008 , 0.62831853])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k])

    def test_trapezoids_and_rectangles(self):
        S_actual = mods.trapezoids_and_rectangles(self.geo)
        S_desired = {
            'pore1': array([1.25821405, 1.37016859]),
            'throat': array([4.e+14, 2.e+00]),
            'pore2': array([1.46368158, 1.53166311])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k])

    def test_pyramids_and_cuboids(self):
        S_actual = mods.pyramids_and_cuboids(self.geo)
        S_desired = {
            'pore1': array([0.82937365, 0.8]),
            'throat': array([1.6e+14, 8.0e-01]),
            'pore2': array([0.85459941, 0.8])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k])

    def test_cubes_and_cuboids(self):
        S_actual = mods.cubes_and_cuboids(self.geo)
        S_desired = {
            'pore1': array([2.4, 1.8]),
            'throat': array([1.6e+14, 8.0e-01]),
            'pore2': array([2.025, 1.4])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k])

    def test_squares_and_rectangles(self):
        S_actual = mods.squares_and_rectangles(self.geo)
        S_desired = {
            'pore1': array([2., 2.]),
            'throat': array([4.e+14, 2.e+00]),
            'pore2': array([2.25, 2.])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k])

    def test_ncylinders_in_series(self):
        S_actual = mods.ncylinders_in_series(self.geo)
        S_desired = {
            'pore1': array([0.49260457, 0.77443401]),
            'throat': array([1.25663706e+14, 4.05813214e-01]),
            'pore2': array([0.563723  , 0.84600371])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k])


if __name__ == '__main__':

    t = ConduitDiffusiveSizeFactorsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
