import openpnm as op
import scipy as sp
import pytest


class BaseTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])

    def teardown_class(self):
        mgr = op.Base.Workspace()
        mgr.clear()


if __name__ == '__main__':

    t = BaseTest()
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
    self = t
