import numpy as np
import scipy as sp
from numpy.testing import assert_allclose

import openpnm as op
import openpnm.models.geometry.throat_surface_area as tsa


class ThroatSurfaceAreaTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.net['throat.diameter'] = np.ones((self.net.Nt, ))
        self.net['throat.length'] = np.ones((self.net.Nt, ))
        self.net['throat.perimeter'] = np.ones((self.net.Nt, ))

    def test_cylinder(self):
        self.net.add_model(propname='throat.surface_area',
                           model=tsa.cylinder, regen_mode="normal")
        assert_allclose(self.net['throat.surface_area'].mean(), np.pi)

    def test_cuboid(self):
        self.net.add_model(propname='throat.surface_area',
                           model=tsa.cuboid, regen_mode="normal")
        assert_allclose(self.net['throat.surface_area'].mean(), 4)

    def test_rectangle(self):
        self.net.add_model(propname='throat.surface_area',
                            model=tsa.rectangle, regen_mode="normal")
        assert_allclose(self.net['throat.surface_area'].mean(), 2)

    def test_extrusion(self):
        self.net.add_model(propname='throat.surface_area',
                           model=tsa.extrusion, regen_mode="normal")
        assert_allclose(self.net['throat.surface_area'].mean(), 1)


if __name__ == '__main__':

    t = ThroatSurfaceAreaTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f"Running test: {item}")
            t.__getattribute__(item)()
