import openpnm as op
import scipy as sp
import pytest


class ProjectTest:

    def setup_class(self):
        self.ws = op.core.Workspace()
        self.ws.clear()

    def test_assign_new_name(self):
        proj = self.ws.new_project()
        new_name = self.ws._gen_name()
        proj.name = new_name
        assert proj.name == new_name
        assert proj.name in self.ws.keys()


if __name__ == '__main__':

    t = ProjectTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
