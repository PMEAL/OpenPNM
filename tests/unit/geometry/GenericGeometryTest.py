import openpnm as op
import pytest


class GenericGeometryTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.geo = op.geometry.SpheresAndCylinders(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)
        self.geo.regenerate_models()

    def teardown_class(self):
        mgr = op.Workspace()
        mgr.clear()

    def test_instantiate_without_network_errors(self):
        with pytest.raises(TypeError):
            op.geometry.SpheresAndCylinders()

    def test_add_locations_after_instantiation(self):
        net = op.network.Cubic(shape=[3, 3, 3])
        geo = op.geometry.SpheresAndCylinders(network=net)
        geo.set_locations(pores=net.Ps)
        assert geo.Np == 27
        assert geo.Nt == 0
        geo.set_locations(throats=net.Ts)
        assert geo.Nt == 54

    def test_plot_histogram(self):
        # Test default args
        self.geo.show_hist()
        # Test with non-default args
        self.geo.show_hist(props=['pore.diameter', 'pore.volume', 'throat.length'])
        # Test layout when num_props = 4 => should be 2 by 2
        self.geo.show_hist(
            props=[
                'pore.diameter',
                'throat.diameter',
                'pore.volume',
                'throat.length'
            ]
        )
        # Test layout when num_props > 4 => should be nrows by 3
        self.geo.show_hist(
            props=[
                'pore.diameter',
                'pore.volume',
                'throat.length',
                'throat.diameter',
                'pore.seed'
            ]
        )


if __name__ == '__main__':

    t = GenericGeometryTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
