import OpenPNM as op
import scipy as sp


class DelaunayVoronoiDualTest:

    def test_default_initialization(self):
        sp.random.seed(seed=0)
        net = op.Network.DelaunayVoronoiDual(num_points=50)
        assert net.Np == 404
        assert net.Nt == 2222
        assert net.num_pores('voronoi') == 281
        assert net.num_pores('delaunay') == 123
        assert net.num_pores('surface') == 114
        assert net.num_pores('boundary') == 73
        assert net.num_throats('voronoi') == 558
        assert net.num_throats('delaunay') == 328
        assert net.num_throats('interconnect') == 1336
        assert net.num_throats('surface') == 185
        assert net.num_throats('boundary') == 443

    def test_cylinder(self):
        sp.random.seed(seed=0)
        net = op.Network.DelaunayVoronoiDual(num_points=50, domain_size=[1, 1])
        assert net.Np == 549
        assert net.Nt == 2963
        assert net.num_pores('voronoi') == 398
        assert net.num_pores('delaunay') == 151
        assert net.num_pores('surface') == 250
        assert net.num_pores('boundary') == 101
        assert net.num_throats('voronoi') == 761
        assert net.num_throats('delaunay') == 414
        assert net.num_throats('interconnect') == 1788
        assert net.num_throats('surface') == 422
        assert net.num_throats('boundary') == 854

    def test_sphere(self):
        sp.random.seed(seed=0)
        net = op.Network.DelaunayVoronoiDual(num_points=50, domain_size=[1])
        assert net.Np == 497
        assert net.Nt == 2806
        assert net.num_pores('voronoi') == 405
        assert net.num_pores('delaunay') == 92
        assert net.num_pores('surface') == 228
        assert net.num_pores('boundary') == 42
        assert net.num_throats('voronoi') == 770
        assert net.num_throats('delaunay') == 416
        assert net.num_throats('interconnect') == 1620
        assert net.num_throats('surface') == 380
        assert net.num_throats('boundary') == 618
