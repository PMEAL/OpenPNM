import openpnm as op


class FickianDiffusionTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[9, 9, 9])
        self.net.add_model_collection(op.models.collections.geometry.spheres_and_cylinders)
        self.net.regenerate_models()
        self.phase = op.phase.Air(network=self.net)
        self.phase.add_model_collection(op.models.collections.physics.basic)
        self.phase.regenerate_models()
        self.phase['throat.diffusive_conductance'] = 1.0

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
