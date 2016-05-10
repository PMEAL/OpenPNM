import scipy as sp
import OpenPNM
import os


class DelaunayTest:

    def setup_class(self):
        bp = sp.array([[0.2, 0.2, 0.2], [0.2, 0.8, 0.2], [0.8, 0.2, 0.2],
                       [0.8, 0.8, 0.2], [0.2, 0.2, 0.8], [0.2, 0.8, 0.8],
                       [0.8, 0.2, 0.8], [0.8, 0.8, 0.8]])
        scale = 1e-4
        self.scale = scale
        p = (sp.random.random([len(bp), 3])-0.5)/10000
        bp += p
        self.mgr = OpenPNM.Base.Workspace()
        self.net = OpenPNM.Network.Delaunay(domain_size=[scale, scale, scale],
                                            base_points=bp*scale)
        self.net.add_boundaries()

    def test_domain_dimensions(self):
        face1 = self.net.pores('bottom_boundary')
        face2 = self.net.pores('top_boundary')
        assert sp.around(self.net.domain_length(face1, face2) / self.scale,
                         3) == 1.0
        assert sp.around(self.net.domain_area(face1) / self.scale**2,
                         3) == 1.0

    def test_export_vor_fibres(self):
        self.net._export_vor_fibres()
        os.remove('fibres.p')
