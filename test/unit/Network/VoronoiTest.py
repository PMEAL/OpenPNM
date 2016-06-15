import OpenPNM as op
import scipy as sp


class VoronoiTest:
    def setup_class(self):
        pass

    def test_reflected(self):
        vor = op.Network.Voronoi(num_cells=50, domain_size=[0.1, 0.1, 0.1],
                                 face_type='reflected')
        assert sp.amax(vor['pore.coords']) < 0.101
        assert sp.amin(vor['pore.coords']) > -0.101

    def test_rough(self):
        vor = op.Network.Voronoi(num_cells=50, domain_size=[0.1, 0.1, 0.1],
                                 face_type='rough')
        assert sp.amax(vor['pore.coords']) < 0.101
        assert sp.amin(vor['pore.coords']) > -0.101
