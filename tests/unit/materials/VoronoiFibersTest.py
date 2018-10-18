import openpnm as op
from openpnm.topotools import reflect_base_points
import scipy as sp
import matplotlib.pyplot as plt


class VoronoiTest:

    def setup_class(self):
        bp = sp.array([[0.25, 0.25, 0.25], [0.25, 0.75, 0.25],
                       [0.75, 0.25, 0.25], [0.75, 0.75, 0.25],
                       [0.75, 0.25, 0.75], [0.25, 0.75, 0.75],
                       [0.25, 0.25, 0.75], [0.75, 0.75, 0.75]])
        scale = 1e-4
        bp = reflect_base_points(bp, [1, 1, 1])*scale
        self.wrk = op.Workspace()
        self.prj = op.materials.VoronoiFibers(fiber_rad=2e-6,
                                              resolution=1e-6,
                                              shape=[scale]*3,
                                              points=bp,
                                              name='test')
        self.net = self.prj.network
        self.del_geom = self.prj.geometries()['test_del']
        self.vor_geom = self.prj.geometries()['test_vor']

    def test_props_all(self):
        a = self.del_geom.props()
        assert sorted(a) == ['pore.area',
                             'pore.centroid',
                             'pore.diameter',
                             'pore.incenter',
                             'pore.indiameter',
                             'pore.vertices',
                             'pore.volume',
                             'throat.area',
                             'throat.c2c',
                             'throat.centroid',
                             'throat.conduit_lengths.pore1',
                             'throat.conduit_lengths.pore2',
                             'throat.conduit_lengths.throat',
                             'throat.diameter',
                             'throat.endpoints.head',
                             'throat.endpoints.tail',
                             'throat.incenter',
                             'throat.indiameter',
                             'throat.length',
                             'throat.normal',
                             'throat.offset_vertices',
                             'throat.perimeter',
                             'throat.shape_factor',
                             'throat.surface_area',
                             'throat.vertices',
                             'throat.volume']

    def test_get_fibre_slice(self):
        slc = self.del_geom._get_fiber_slice(index=[0, 50, 0])
        assert sp.shape(slc) == (101, 101)

    def test_plot_pore(self):
        self.del_geom.plot_pore(pores=self.del_geom.pores())
        plt.close('all')

    def test_plot_throat(self):
        self.del_geom.plot_throat(throats=[0])
        plt.close('all')

    def test_vertex_dimension(self):
        prj = op.materials.VoronoiFibers(num_points=50,
                                         fiber_rad=0.2,
                                         resolution=0.1,
                                         shape=[3, 2, 1],
                                         name='test2')
        net = prj.network
        del_geom = prj.geometries()['test2_del']
        B1 = net.pores(['left', 'delaunay'], mode='xnor')
        B2 = net.pores(['right', 'delaunay'], mode='xnor')
        assert del_geom.vertex_dimension(B1, B2, 'volume') == 6.0
        assert del_geom.vertex_dimension(B1, B2, 'area') == 3.0
        assert del_geom.vertex_dimension(B1, B2, 'length') == 2.0
        assert del_geom.vertex_dimension(B1, B2, 'area_xy') == 6.0
        assert del_geom.vertex_dimension(B1, B2, 'area_yz') == 2.0
        assert del_geom.vertex_dimension(B1, B2, 'area_xz') == 3.0
        assert del_geom.vertex_dimension(B1, B2, 'minmax') == \
            [0.0, 3.0, 0.0, 2.0, 0.0, 1.0]

    def test_linear_scale(self):
        prj = op.materials.VoronoiFibers(num_points=50,
                                         fiber_rad=0.2,
                                         resolution=0.1,
                                         shape=[2, 2, 2],
                                         name='test3',
                                         linear_scale=[1, 1, 2])
        net = prj.network
        del_geom = prj.geometries()['test3_del']
        B1 = net.pores(['top', 'delaunay'], mode='xnor')
        B2 = net.pores(['bottom', 'delaunay'], mode='xnor')
        assert del_geom.vertex_dimension(B1, B2, 'volume') == 8.0
        assert del_geom.vertex_dimension(B1, B2, 'length') == 2.0

    def test_linear_scale_wrong_shape(self):
        prj = op.materials.VoronoiFibers(num_points=50,
                                         fiber_rad=0.2,
                                         resolution=0.1,
                                         shape=[2, 2, 2],
                                         name='test4',
                                         linear_scale=[1, 1])
        net = prj.network
        del_geom = prj.geometries()['test4_del']
        B1 = net.pores(['top', 'delaunay'], mode='xnor')
        B2 = net.pores(['bottom', 'delaunay'], mode='xnor')
        assert del_geom.vertex_dimension(B1, B2, 'volume') == 8.0
        assert del_geom.vertex_dimension(B1, B2, 'length') == 2.0


if __name__ == '__main__':

    t = VoronoiTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t
