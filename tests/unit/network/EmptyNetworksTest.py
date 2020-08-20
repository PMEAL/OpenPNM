import openpnm as op


class EmptyNetworksTest:
    def setup_class(self):
        pass

    def teardown_class(self):
        pass

    def test_each(self):
        pn = op.network.Cubic()
        pn = op.network.Bravais()
        pn = op.network.CubicDual()
        pn = op.network.CubicTemplate()
        pn = op.network.Delaunay()
        pn = op.network.DelaunayVoronoiDual()
        pn = op.network.Gabriel()
        pn = op.network.Voronoi()


if __name__ == '__main__':

    t = EmptyNetworksTest()
    t.setup_class()
    self = t
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
