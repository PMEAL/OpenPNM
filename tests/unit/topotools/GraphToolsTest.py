import openpnm as op
import numpy as np
from numpy.testing import assert_approx_equal
from openpnm import topotools


class GraphToolsTest:

    def setup_class(self):
        self.ws = op.Workspace()

    def teardown_class(self):
        self.ws.clear()


if __name__ == '__main__':

    t = GraphToolsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
