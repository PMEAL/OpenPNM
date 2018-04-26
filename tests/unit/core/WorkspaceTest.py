import openpnm as op
import scipy as sp
import pytest


class WorkspaceTest:

    def setup_class(self):
        self.ws = op.core.Workspace()
        self.ws.clear()
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.phase = op.phases.Air(network=self.net)

    def test_new_project_no_name(self):
        proj = self.ws.new_project()
        assert proj.name in self.ws.keys()

    def test_new_project_w_name(self):
        proj = self.ws.new_project(name='bob')
        assert proj.name in self.ws.keys()

    def test_new_project_duplicate_name(self):
        proj = self.ws.new_project('foo')
        with pytest.raises(Exception):
            proj = self.ws.new_project('foo')
        assert proj.name in self.ws.keys()

    def test_assign_project(self):
        proj = self.ws.new_project()
        with pytest.raises(Exception):
            self.ws[proj.name] = proj
        old_name = proj.name
        new_name = self.ws._gen_name()
        self.ws[new_name] = proj
        assert proj.name == new_name
        assert proj.name in self.ws.keys()
        assert old_name not in self.ws.keys()

    def test_str(self):
        s = self.ws.__str__().split('\n')
        assert 'OpenPNM Version' in s[1]


if __name__ == '__main__':

    t = WorkspaceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
