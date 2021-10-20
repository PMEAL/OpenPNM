import openpnm as op


class FickianDiffusionTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[9, 9, 9])
        self.geo = op.geometry._StickAndBall(network=self.net,
                                             pores=self.net.Ps,
                                             throats=self.net.Ts)
        self.phase = op.phases.Air(network=self.net)
        self.phase['pore.mole_fraction'] = 0
        self.phys = op.physics.GenericPhysics(network=self.net,
                                              phase=self.phase,
                                              geometry=self.geo)
        self.phys['throat.diffusive_conductance'] = 1.0

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()


if __name__ == '__main__':

    t = FickianDiffusionTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t
