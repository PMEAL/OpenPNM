import openpnm as op
import pytest


class PhaseTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[10, 10, 10])

    def teardown_class(self):
        mgr = op.Workspace()
        mgr.clear()

    def test_instantiate_without_network_fails(self):
        with pytest.raises(TypeError):
            op.phase.Phase()


if __name__ == '__main__':

    t = PhaseTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
