import OpenPNM
import scipy as sp
import matplotlib.pyplot as plt
from OpenPNM.Geometry import models as gm
import OpenPNM.Utilities.misc as misc


class VoronoiTest:

    def setup_class(self):
        bp = sp.array([[0.2, 0.2, 0.2], [0.2, 0.8, 0.2], [0.8, 0.2, 0.2],
                       [0.8, 0.8, 0.2], [0.2, 0.2, 0.8], [0.2, 0.8, 0.8],
                       [0.8, 0.2, 0.8], [0.8, 0.8, 0.8]])
        scale = 1e-4
        p = (sp.random.random([len(bp), 3])-0.5)/10000
        bp += p
        self.ctrl = OpenPNM.Base.Controller()
        self.net = OpenPNM.Network.Delaunay(domain_size=[scale, scale, scale],
                                            base_points=bp*scale)
        self.net.add_boundaries()
        Ps = self.net.pores()
        Ts = self.net.throats()
        fibre_rad = 5e-6
        self.geo_vox = OpenPNM.Geometry.Voronoi(network=self.net,
                                                pores=Ps,
                                                throats=Ts,
                                                fibre_rad=fibre_rad,
                                                voxel_vol=True,
                                                name='vor_vox')

    def test_props_all(self):
        a = self.geo_vox.props()
        assert sorted(a) == ['pore.area', 'pore.centroid', 'pore.diameter',
                             'pore.fibre_voxels', 'pore.indiameter',
                             'pore.pore_voxels', 'pore.seed',
                             'pore.vertices', 'pore.volume',
                             'throat.area', 'throat.c2c', 'throat.centroid',
                             'throat.diameter', 'throat.incentre',
                             'throat.indiameter', 'throat.length',
                             'throat.normal', 'throat.offset_vertices',
                             'throat.perimeter', 'throat.seed',
                             'throat.shape_factor', 'throat.surface_area',
                             'throat.vertices', 'throat.volume']

    def test_get_fibre_slice(self):
        slc = self.geo_vox.get_fibre_slice(index=[0, 50, 0])
        assert sp.shape(slc) == (101, 101)

    def test_compress_geom(self):
        b1 = self.net.pores('bottom_boundary')
        b2 = self.net.pores('top_boundary')
        height1 = self.net.domain_length(b1, b2)
        self.geo_vox.compress_geometry(factor=[1, 1, 0.5])
        height2 = self.net.domain_length(b1, b2)
        assert sp.around(height1/height2, 5) == 2.0
    
    def test_plot_porosity_profile(self):
        fig = plt.figure()
        self.geo_vox.plot_porosity_profile(fig)
        del fig
    
    def test_plot_fibre_slice(self):
        fig = plt.figure()
        self.geo_vox.plot_fibre_slice(fig)
        del fig
    
    def test_export_fibre_image(self):
        self.geo_vox.export_fibre_image(mat_file='OpenPNMFibres')
    
    def test_make_fibre_image(self):
        del(self.geo_vox._fibre_image)
        self.geo_vox.make_fibre_image()
        assert hasattr(self.geo_vox,'_fibre_image') == True
    
    def test_voronoi_vert(self):
        self.ctrl.purge_object(self.geo_vox)
        fibre_rad = 5e-6
        Ps = self.net.pores()
        Ts = self.net.throats()
        self.geo_vert = OpenPNM.Geometry.Voronoi(network=self.net,
                                                pores=Ps,
                                                throats=Ts,
                                                fibre_rad=fibre_rad,
                                                voxel_vol=False,
                                                name='vor_vert')
        self.geo_vert.models.add(propname = 'throat.area2',
                                 model=gm.throat_area.voronoi)
        assert sp.all(sp.absolute(1 - self.geo_vert["throat.area2"]/
                                  self.geo_vert['throat.area']) < 0.05)
        self.geo_vert.models.add(propname = 'throat.centroid2',
                                 model=gm.throat_centroid.centre_of_mass)
        assert sp.absolute(sp.sum(self.geo_vert['throat.centroid']-
                                  self.geo_vert['throat.centroid2'])) < 1e-10
        self.geo_vert.models.add(propname = 'throat.perimeter2',
                                 model=gm.throat_perimeter.voronoi)
        assert sp.all(sp.absolute(1 - self.geo_vert['throat.perimeter2']/
                                  self.geo_vert['throat.perimeter']) < 0.05)
        self.geo_vert.models.add(propname = 'throat.length2',
                                 model=gm.throat_length.voronoi)
        cond_lengths = sp.sum(misc.conduit_lengths(network = self.net,
                                                   mode='centroid'),axis=1)
        assert sp.all(sp.absolute(1 - cond_lengths/
                      self.geo_vert['throat.length2']) < 0.05)
