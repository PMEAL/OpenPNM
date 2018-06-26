import openpnm as op
from openpnm.topotools import reflect_base_points
import openpnm.models.geometry as gm
import scipy as sp


class VoronoiTest:

    def setup_class(self):
        bp = sp.array([[0.25, 0.25, 0.25], [0.25, 0.75, 0.25],
                       [0.75, 0.25, 0.25], [0.75, 0.75, 0.25],
                       [0.75, 0.25, 0.75], [0.25, 0.75, 0.75],
                       [0.25, 0.25, 0.75], [0.75, 0.75, 0.75]])
        scale = 1e-4
        bp = reflect_base_points(bp, [1, 1, 1])*scale
        self.wrk = op.Workspace()
        self.net = op.materials.VoronoiFibers(fiber_rad=2e-6,
                                              resolution=1e-6,
                                              shape=[scale]*3,
                                              points=bp,
                                              name='test')
        self.prj = self.net.project
        self.del_geom = self.prj.geometries()['test_del']
        self.vor_geom = self.prj.geometries()['test_vor']

    def test_props_all(self):
        a = self.del_geom.props()
        assert sorted(a) == ['pore.area', 'pore.centroid', 'pore.diameter',
                             'pore.incenter', 'pore.indiameter',
                             'pore.vertices', 'pore.volume',
                             'throat.area', 'throat.c2c', 'throat.centroid',
                             'throat.diameter', 'throat.incenter',
                             'throat.indiameter', 'throat.length',
                             'throat.normal', 'throat.offset_vertices',
                             'throat.perimeter',
                             'throat.shape_factor', 'throat.surface_area',
                             'throat.vertices', 'throat.volume']

    def test_get_fibre_slice(self):
        slc = self.del_geom._get_fiber_slice(index=[0, 50, 0])
        assert sp.shape(slc) == (101, 101)


if __name__ == '__main__':

    t = VoronoiTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t
