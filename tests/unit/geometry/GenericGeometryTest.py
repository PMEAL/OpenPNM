import openpnm as op


class GenericGeometryTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.geo = op.geometry.StickAndBall(network=self.net,
                                            pores=self.net.Ps,
                                            throats=self.net.Ts)
        self.geo.regenerate_models()

    def teardown_class(self):
        mgr = op.Workspace()
        mgr.clear()

    def test_plot_histogram(self):
        self.geo.show_hist()
        self.geo.show_hist(props=['pore.diameter', 'pore.volume',
                                  'throat.length'])
        self.geo.show_hist(props=['pore.diameter', 'pore.volume',
                                  'throat.length', 'throat.diameter',
                                  'pore.seed'])


if __name__ == '__main__':

    t = GenericGeometryTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
