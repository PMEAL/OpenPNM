from numpy import array
from numpy.testing import assert_allclose
import openpnm as op
import numpy as np
import openpnm.models.geometry.diffusive_size_factors as mods
import openpnm.models.geometry as gm


class ConduitDiffusiveSizeFactorsTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 1, 1], spacing=1.0)
        self.net.regenerate_models()
        self.net["throat.diameter"] = 0.4
        self.net["pore.diameter"] = [1.2, 0.9, 0.7]
        self.air = op.phase.Air(network=self.net)
        self.air.regenerate_models()

    def test_spheres_and_cylinders(self):
        SF = mods.spheres_and_cylinders(self.net)
        S_actual = {'pore1': SF[:, 0],
                    'throat': SF[:, 1],
                    'pore2': SF[:, 2]}
        S_desired = {
            'pore1': array([0.93875728, 0.97459088]),
            'throat': array([1.25663706e+14, 4.05813214e-01]),
            'pore2': array([0.82884551, 0.94886745])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k])

    def test_circles_and_rectangles(self):
        SF = mods.circles_and_rectangles(self.net)
        S_actual = {'pore1': SF[:, 0],
                    'throat': SF[:, 1],
                    'pore2': SF[:, 2]}
        S_desired = {
            'pore1': array([1.53390798, 1.80140852]),
            'throat': array([4.00000000e+14, 1.29174358e+00]),
            'pore2': array([1.65097536, 2.07781253])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k])

    def test_cones_and_cylinders(self):
        SF = mods.cones_and_cylinders(self.net)
        S_actual = {'pore1': SF[:, 0],
                    'throat': SF[:, 1],
                    'pore2': SF[:, 2]}
        S_desired = {
            'pore1': array([0.65138854, 0.62831853]),
            'throat': array([1.25663706e+14, 6.28318531e-01]),
            'pore2': array([0.6712008, 0.62831853])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k])

    def test_trapezoids_and_rectangles(self):
        SF = mods.trapezoids_and_rectangles(self.net)
        S_actual = {'pore1': SF[:, 0],
                    'throat': SF[:, 1],
                    'pore2': SF[:, 2]}
        S_desired = {
            'pore1': array([1.25821405, 1.37016859]),
            'throat': array([4.e+14, 2.e+00]),
            'pore2': array([1.46368158, 1.53166311])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k])

    def test_pyramids_and_cuboids(self):
        SF = mods.pyramids_and_cuboids(self.net)
        S_actual = {'pore1': SF[:, 0],
                    'throat': SF[:, 1],
                    'pore2': SF[:, 2]}
        S_desired = {
            'pore1': array([0.82937365, 0.8]),
            'throat': array([1.6e+14, 8.0e-01]),
            'pore2': array([0.85459941, 0.8])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k])

    def test_cubes_and_cuboids(self):
        SF = mods.cubes_and_cuboids(self.net)
        S_actual = {'pore1': SF[:, 0],
                    'throat': SF[:, 1],
                    'pore2': SF[:, 2]}
        S_desired = {
            'pore1': array([2.4, 1.8]),
            'throat': array([1.6e+14, 8.0e-01]),
            'pore2': array([2.025, 1.4])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k])

    def test_squares_and_rectangles(self):
        SF = mods.squares_and_rectangles(self.net)
        S_actual = {'pore1': SF[:, 0],
                    'throat': SF[:, 1],
                    'pore2': SF[:, 2]}
        S_desired = {
            'pore1': array([2., 2.]),
            'throat': array([4.e+14, 2.e+00]),
            'pore2': array([2.25, 2.])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k])

    def test_ncylinders_in_series(self):
        SF = mods.ncylinders_in_series(self.net)
        S_actual = {'pore1': SF[:, 0],
                    'throat': SF[:, 1],
                    'pore2': SF[:, 2]}
        S_desired = {
            'pore1': array([0.49260457, 0.77443401]),
            'throat': array([1.25663706e+14, 4.05813214e-01]),
            'pore2': array([0.563723, 0.84600371])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k])

    def test_intersecting_cones(self):
        self.net.add_model(propname='throat.coords',
                           model=gm.throat_centroid.pore_coords)
        SF = mods.intersecting_cones(self.net)
        S_actual = {'pore1': SF[:, 0],
                    'throat': SF[:, 1],
                    'pore2': SF[:, 2]}
        S_desired = {
            'pore1': array([0.75398224, 0.56548668]),
            'throat': array([np.inf, np.inf]),
            'pore2': array([0.56548668, 0.43982297])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k])

    def test_hybrid_cones_and_cylinders(self):
        self.net.add_model(propname='throat.coords',
                           model=gm.throat_centroid.pore_coords)
        self.net["pore.diameter"][0] = 1.6
        SF = mods.hybrid_cones_and_cylinders(self.net)
        S_actual = {'pore1': SF[:, 0],
                    'throat': SF[:, 1],
                    'pore2': SF[:, 2]}
        S_desired = {
            'pore1': array([1.00530965, 0.62831853]),
            'throat': array([np.inf, 0.62831853]),
            'pore2': array([0.56548668, 0.62831853])
        }
        for k, v in S_actual.items():
            assert_allclose(v, S_desired[k])
        self.net["pore.diameter"][0] = 1.2


if __name__ == '__main__':

    t = ConduitDiffusiveSizeFactorsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
