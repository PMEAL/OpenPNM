import OpenPNM
import scipy as sp
import OpenPNM.Utilities.vertexops as vo


class VertexOpsTest:
    def setup_class(self):
        bp = sp.array([[0.2, 0.2, 0.2], [0.2, 0.8, 0.2], [0.8, 0.2, 0.2],
                       [0.8, 0.8, 0.2], [0.2, 0.2, 0.8], [0.2, 0.8, 0.8],
                       [0.8, 0.2, 0.8], [0.8, 0.8, 0.8]])
        scale = 1e-4
        p = (sp.random.random([len(bp), 3])-0.5)/10000
        bp += p
        self.mgr = OpenPNM.Base.Workspace()
        self.net = OpenPNM.Network.Delaunay(domain_size=[scale, scale, scale],
                                            base_points=bp*scale)
        self.net.add_boundaries()
        Ps = self.net.pores()
        Ts = self.net.throats()
        self.fibre_rad = 5e-6
        self.geo = OpenPNM.Geometry.Voronoi(network=self.net,
                                            pores=Ps,
                                            throats=Ts,
                                            fibre_rad=self.fibre_rad,
                                            voxel_vol=False,
                                            name='vor')

    def test_scale(self):
        factor = [1, 1, 0.5]
        vo.scale(network=self.net,
                 scale_factor=factor,
                 linear_scaling=[True, False, False],
                 preserve_vol=False)

    def test_porosity(self):
        por = vo.porosity(self.net)
        assert por < 1.0

    def test_tortuosity(self):
        tor = vo.tortuosity(self.net)
        assert sp.all(tor > 1.0)

    def test_pore2centroid(self):
        temp_coords = self.net['pore.coords']
        self.geo['pore.centroid'] = sp.ones([self.geo.num_pores(), 3])
        vo.pore2centroid(self.net)
        assert sp.sum(self.net['pore.coords'] -
                      sp.ones([self.geo.num_pores(), 3])) == 0.0
        self.net['pore.coords'] = temp_coords

    # def test_plot_pore(self):
    #    vo.plot_pore(self.geo, self.geo.pores())

    # def test_plot_throat(self):
    #    vo.plot_throat(self.geo, [1, 2, 3, 4])

    def test_rotate_and_chop(self):
        throat_verts = self.geo["throat.vertices"][0]
        throat_normal = self.geo["throat.normal"][0]
        test = vo.rotate_and_chop(throat_verts, throat_normal, [0, 1, 0])
        r, c = sp.shape(test)
        assert r == len(throat_verts)
        assert c == 2
