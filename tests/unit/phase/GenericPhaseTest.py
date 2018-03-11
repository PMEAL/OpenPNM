import openpnm as op
import scipy as sp


class GenericPhaseTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[10, 10, 10])

    def teardown_class(self):
        mgr = op.core.Workspace()
        mgr.clear()
