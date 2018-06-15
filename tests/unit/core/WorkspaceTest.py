import openpnm as op
import scipy as sp
import pytest
import os
import pickle


class WorkspaceTest:

    def setup_class(self):
        self.ws = op.core.Workspace()
        self.ws.clear()
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.phase = op.phases.Air(network=self.net)

    def test_new_project_no_name(self):
        proj = self.ws.new_project()
        assert proj.name in self.ws.keys()
        self.ws.clear()

    def test_new_project_w_name(self):
        proj = self.ws.new_project(name='bob')
        assert proj.name in self.ws.keys()
        self.ws.clear()

    def test_new_project_duplicate_name(self):
        proj = self.ws.new_project('foo')
        with pytest.raises(Exception):
            proj = self.ws.new_project('foo')
        assert proj.name in self.ws.keys()
        self.ws.clear()

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
        self.ws.clear()

    def test_str(self):
        proj = self.ws.new_project()
        op.network.Cubic(shape=[3, 3, 3], project=proj)
        s = self.ws.__str__().split('\n')
        assert 'OpenPNM Version' in s[1]
        self.ws.clear()

    def test_save_and_load_project(self):
        proj = self.ws.new_project('test_proj')
        net = op.network.Cubic(shape=[3, 3, 3], project=proj)
        op.phases.Air(network=net)
        self.ws.save_project(proj)
        assert proj.name in self.ws.keys()
        self.ws.close_project(proj)
        assert 'test_proj' not in self.ws.keys()
        assert proj.workspace == {}
        self.ws.load_project(filename='test_proj.pnm')
        assert 'test_proj' in self.ws.keys()
        self.ws.clear()
        os.remove('test_proj.pnm')

    def test_save_and_load_project_from_pickled_list(self):
        proj = self.ws.new_project()
        pn = op.network.Cubic(shape=[3, 3, 3], project=proj)
        air = op.phases.Air(network=pn)
        pickle.dump([pn, air], open('test.pnm', 'wb'))
        self.ws.clear()
        self.ws.load_project('test.pnm')
        self.ws.clear()
        os.remove('test.pnm')

    def test_save_and_load_project_from_pickled_object(self):
        a = sp.ones((10, ))
        pickle.dump(a, open('single_object.pnm', 'wb'))
        self.ws.clear()
        with pytest.warns(UserWarning):
            self.ws.load_project('single_object.pnm')
        b = {'test': a}
        pickle.dump(b, open('single_object.pnm', 'wb'))
        self.ws.clear()
        with pytest.warns(UserWarning):
            self.ws.load_project('single_object.pnm')
        os.remove('single_object.pnm')

    def test_load_project_with_name_conflict(self):
        self.ws.clear()
        proj = self.ws.new_project(name='test')
        pn = op.network.Cubic(shape=[3, 3, 3], project=proj)
        op.phases.Air(network=pn)
        self.ws.save_project(proj, filename='test.pnm')
        with pytest.warns(UserWarning):
            self.ws.load_project('test.pnm')
        os.remove('test.pnm')

    def test_save_and_load_workspace(self):
        proj1 = self.ws.new_project('test_proj_1')
        proj2 = self.ws.new_project('test_proj_2')
        op.network.Cubic(shape=[3, 3, 3], project=proj1, name='net1')
        op.network.Cubic(shape=[3, 3, 3], project=proj2, name='net2')
        self.ws.save_workspace(filename='workspace_test')
        self.ws.clear()
        assert 'test_proj_1' not in self.ws.keys()
        assert 'test_proj_2' not in self.ws.keys()
        self.ws.load_workspace('workspace_test')
        assert 'test_proj_1' in self.ws.keys()
        assert 'test_proj_2' in self.ws.keys()
        self.ws.clear()
        os.remove('workspace_test.pnm')


if __name__ == '__main__':

    t = WorkspaceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
