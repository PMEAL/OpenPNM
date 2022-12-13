import numpy as np
from numpy.testing import assert_array_almost_equal

import openpnm as op
import openpnm.models.geometry as gm
import openpnm.models.misc as mm


class ThroatLengthTest:

    def setup_class(self):
        np.random.seed(0)
        self.net = op.network.Demo(shape=[5, 1, 1], spacing=1)
        self.net.add_model(propname='pore.seed',
                           model=mm.constant,
                           value=0.5)
        self.net.regenerate_models()
        self.net['throat.coords'] = np.mean(self.net.coords[self.net.conns], axis=1)

    def test_circles_and_rectangles(self):
        del self.net['throat.length']
        self.net.add_model(propname='throat.length',
                           model=gm.throat_length.circles_and_rectangles)
        assert_array_almost_equal(self.net['throat.length'], 0.5669873)

    def test_trapezoids_and_rectangles(self):
        del self.net['throat.length']
        self.net.add_model(propname='throat.length',
                           model=gm.throat_length.trapezoids_and_rectangles)
        assert_array_almost_equal(self.net['throat.length'], 0.5)

    def test_cubes_and_cuboids(self):
        del self.net['throat.length']
        self.net.add_model(propname='throat.length',
                           model=gm.throat_length.cubes_and_cuboids)
        assert_array_almost_equal(self.net['throat.length'], 0.5)
        # Make sure it still works if throat.spacing is not defined
        del self.net['throat.spacing']
        self.net.add_model(propname='throat.length',
                           model=gm.throat_length.cubes_and_cuboids)
        assert_array_almost_equal(self.net['throat.length'], 0.5)
        # Add back throat.spacing so other tests don't fail
        self.net.regenerate_models()

    def test_squares_and_rectangles(self):
        del self.net['throat.length']
        self.net.add_model(propname='throat.length',
                           model=gm.throat_length.squares_and_rectangles)
        assert_array_almost_equal(self.net['throat.length'], 0.5)

    def test_intersecting_cones(self):
        del self.net['throat.length']
        self.net.add_model(propname='throat.length',
                           model=gm.throat_length.intersecting_cones)
        assert_array_almost_equal(self.net['throat.length'], 0.0)

    def test_intersecting_pyramids(self):
        del self.net['throat.length']
        self.net.add_model(propname='throat.length',
                           model=gm.throat_length.intersecting_pyramids)
        assert_array_almost_equal(self.net['throat.length'], 0.0)

    def test_hybrid_cones_and_cylinders(self):
        del self.net['throat.length']
        self.net.add_model(propname='throat.length',
                           model=gm.throat_length.cones_and_cylinders)
        assert_array_almost_equal(self.net['throat.length'], 0.5)
        # Make sure it still works if throat.spacing is not defined
        del self.net['throat.spacing']
        self.net.add_model(propname='throat.length',
                           model=gm.throat_length.cones_and_cylinders)
        assert_array_almost_equal(self.net['throat.length'], 0.5)
        # Add back throat.spacing so other tests don't fail
        self.net.regenerate_models()

    def test_hybrid_trapezoids_and_rectangles(self):
        del self.net['throat.length']
        self.net.add_model(propname='throat.length',
                           model=gm.throat_length.hybrid_trapezoids_and_rectangles)
        assert_array_almost_equal(self.net['throat.length'], 0.5)

    def test_hybrid_pyramids_and_cuboids(self):
        del self.net['throat.length']
        self.net.add_model(propname='throat.length',
                           model=gm.throat_length.hybrid_pyramids_and_cuboids)
        assert_array_almost_equal(self.net['throat.length'], 0.5)



if __name__ == '__main__':

    t = ThroatLengthTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
