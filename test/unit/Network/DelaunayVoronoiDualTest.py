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
        assert net.num_pores('surface') == 190
        assert net.num_pores('internal') == 215
        assert net.num_throats('voronoi') == 556
        assert net.num_throats('delaunay') == 327
        assert net.num_throats('interconnect') == 1334
        assert net.num_throats('surface') == 564
        assert net.num_throats('internal') == 1653
        Np_int = net.num_pores('internal')
        Np_surf = net.num_pores('surface')
        assert (Np_int + Np_surf) == net.Np
        Nt_int = net.num_throats('internal')
        Nt_surf = net.num_throats('surface')
        assert (Nt_int + Nt_surf) == net.Nt
        Np_vor = net.num_pores(labels=['voronoi', 'surface'],
                               mode='intersection')
        Np_del = net.num_pores(labels=['delaunay', 'surface'],
                               mode='intersection')
        assert (Np_vor + Np_del) == net.num_pores('surface')
        Nt_vor = net.num_throats(labels=['voronoi', 'surface'],
                                 mode='intersection')
        Nt_del = net.num_throats(labels=['delaunay', 'surface'],
                                 mode='intersection')
        Nt_int = net.num_throats(labels=['interconnect', 'surface'],
                                 mode='intersection')
        assert (Nt_vor + Nt_del + Nt_int) == net.num_throats('surface')

    def test_cylinder(self):
        sp.random.seed(seed=0)
        net = op.Network.DelaunayVoronoiDual(num_points=50, domain_size=[1, 1])
        assert net.Np == 540
        assert net.Nt == 2882
        assert net.num_pores('voronoi') == 388
        assert net.num_pores('delaunay') == 152
        assert net.num_pores('surface') == 240
        assert net.num_pores('internal') == 300
        assert net.num_throats('voronoi') == 739
        assert net.num_throats('delaunay') == 401
        assert net.num_throats('interconnect') == 1742
        assert net.num_throats('surface') == 663
        assert net.num_throats('internal') == 2219
        Np_int = net.num_pores('internal')
        Np_surf = net.num_pores('surface')
        assert (Np_int + Np_surf) == net.Np
        Nt_int = net.num_throats('internal')
        Nt_surf = net.num_throats('surface')
        assert (Nt_int + Nt_surf) == net.Nt
        Np_vor = net.num_pores(labels=['voronoi', 'surface'],
                               mode='intersection')
        Np_del = net.num_pores(labels=['delaunay', 'surface'],
                               mode='intersection')
        assert (Np_vor + Np_del) == net.num_pores('surface')
        Nt_vor = net.num_throats(labels=['voronoi', 'surface'],
                                 mode='intersection')
        Nt_del = net.num_throats(labels=['delaunay', 'surface'],
                                 mode='intersection')
        Nt_int = net.num_throats(labels=['interconnect', 'surface'],
                                 mode='intersection')
        assert (Nt_vor + Nt_del + Nt_int) == net.num_throats('surface')

    def test_sphere(self):
        sp.random.seed(seed=0)
        net = op.Network.DelaunayVoronoiDual(num_points=50, domain_size=[1])
        assert net.Np == 490
        assert net.Nt == 2766
        assert net.num_pores('voronoi') == 399
        assert net.num_pores('delaunay') == 91
        assert net.num_pores('surface') == 109
        assert net.num_pores('internal') == 381
        assert net.num_throats('voronoi') == 759
        assert net.num_throats('delaunay') == 411
        assert net.num_throats('interconnect') == 1596
        assert net.num_throats('surface') == 255
        assert net.num_throats('internal') == 2511
        Np_int = net.num_pores('internal')
        Np_surf = net.num_pores('surface')
        assert (Np_int + Np_surf) == net.Np
        Nt_int = net.num_throats('internal')
        Nt_surf = net.num_throats('surface')
        assert (Nt_int + Nt_surf) == net.Nt
        Np_vor = net.num_pores(labels=['voronoi', 'surface'],
                               mode='intersection')
        Np_del = net.num_pores(labels=['delaunay', 'surface'],
                               mode='intersection')
        assert (Np_vor + Np_del) == net.num_pores('surface')
        Nt_vor = net.num_throats(labels=['voronoi', 'surface'],
                                 mode='intersection')
        Nt_del = net.num_throats(labels=['delaunay', 'surface'],
                                 mode='intersection')
        Nt_int = net.num_throats(labels=['interconnect', 'surface'],
                                 mode='intersection')
        assert (Nt_vor + Nt_del + Nt_int) == net.num_throats('surface')
