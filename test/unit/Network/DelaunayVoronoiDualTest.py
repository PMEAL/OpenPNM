import OpenPNM as op
import scipy as sp


class DelaunayVoronoiDualTest:

    def test_default_initialization(self):
        sp.random.seed(seed=0)
        net = op.Network.DelaunayVoronoiDual(num_points=50)
        assert net.Np == 405
        assert net.Nt == 2217
        assert net.num_pores('voronoi') == 280
        assert net.num_pores('delaunay') == 125
        assert net.num_pores('surface') == 115
        assert net.num_pores('boundary') == 75
        assert net.num_throats('voronoi') == 556
        assert net.num_throats('delaunay') == 327
        assert net.num_throats('interconnect') == 1334
        assert net.num_throats('surface') == 188
        assert net.num_throats('boundary') == 451

    def test_cylinder(self):
        sp.random.seed(seed=0)
        net = op.Network.DelaunayVoronoiDual(num_points=50, domain_size=[1, 1])
        assert net.Np == 540
        assert net.Nt == 2882
        assert net.num_pores('voronoi') == 388
        assert net.num_pores('delaunay') == 152
        assert net.num_pores('surface') == 229
        assert net.num_pores('boundary') == 102
        assert net.num_throats('voronoi') == 739
        assert net.num_throats('delaunay') == 401
        assert net.num_throats('interconnect') == 1742
        assert net.num_throats('surface') == 384
        assert net.num_throats('boundary') == 803

    def test_sphere(self):
        sp.random.seed(seed=0)
        net = op.Network.DelaunayVoronoiDual(num_points=50, domain_size=[1])
        assert net.Np == 490
        assert net.Nt == 2766
        assert net.num_pores('voronoi') == 399
        assert net.num_pores('delaunay') == 91
        assert net.num_pores('surface') == 232
        assert net.num_pores('boundary') == 41
        assert net.num_throats('voronoi') == 759
        assert net.num_throats('delaunay') == 411
        assert net.num_throats('interconnect') == 1596
        assert net.num_throats('surface') == 388
        assert net.num_throats('boundary') == 626
