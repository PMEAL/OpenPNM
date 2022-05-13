import os
import pytest
import openpnm as op
from numpy.testing import assert_allclose


class WorkspaceTest:

    def setup_class(self):
        self.ws = op.Workspace()
        self.ws.clear()
        self.net = op.network.Cubic(shape=[5, 5, 5])
        self.phase = op.phase.Air(network=self.net)

    def test_new_project_no_name(self):
        proj = self.ws.new_project()
        assert proj.name in self.ws.keys()
        self.ws.clear()

    def test_new_project_w_name(self):
        proj = self.ws.new_project(name='bob')
        assert proj.name in self.ws.keys()
        self.ws.clear()

    def test_new_project_duplicate_name(self):
        proj1 = self.ws.new_project('foo')
        proj2 = self.ws.new_project('foo')
        assert proj1.name != proj2.name
        assert proj2.name in self.ws.keys()
        self.ws.clear()

    # def test_assign_project(self):
    #     proj = self.ws.new_project()
    #     with pytest.raises(Exception):
    #         self.ws.new_project(name=proj.name)
    #     old_name = proj.name
    #     new_name = self.ws._validate_name()
    #     self.ws[new_name] = proj
    #     assert proj.name == new_name
    #     assert proj.name in self.ws.keys()
    #     assert old_name not in self.ws.keys()
    #     self.ws.clear()

    def test_str(self):
        proj = self.ws.new_project()
        op.network.Cubic(shape=[3, 3, 3], project=proj)
        s = self.ws.__str__().split('\n')
        assert len(s) == 8
        self.ws.clear()

    # def test_save_and_load_project(self):
    #     proj = self.ws.new_project('test_proj')
    #     net = op.network.Cubic(shape=[3, 3, 3], project=proj)
    #     op.phase.Air(network=net)
    #     self.ws.save_project(proj)
    #     assert proj.name in self.ws.keys()
    #     self.ws.close_project(proj)
    #     assert 'test_proj' not in self.ws.keys()
    #     assert proj.workspace == {}
    #     proj = self.ws.load_project(filename='test_proj.pnm')
    #     assert 'test_proj' in self.ws.keys()
    #     assert isinstance(proj, op.Project)
    #     shape = op.topotools.get_shape(proj.network)
    #     assert_allclose(shape, [3, 3, 3])
    #     self.ws.clear()
    #     try:
    #         os.remove('test_proj.pnm')
    #     except PermissionError:
    #         print('Could not delete test_proj.pnm')

    # def test_load_project_with_name_conflict(self):
    #     self.ws.clear()
    #     proj = self.ws.new_project(name='test')
    #     pn = op.network.Cubic(shape=[3, 3, 3], project=proj)
    #     op.phase.Air(network=pn)
    #     self.ws.save_project(proj, filename='test.pnm')
    #     self.ws.load_project('test.pnm')
    #     assert set(self.ws.keys()) == set(['test', 'proj_01'])
    #     os.remove('test.pnm')

    # def test_save_and_load_project_from_pickled_list(self):
    #     proj = self.ws.new_project()
    #     pn = op.network.Cubic(shape=[3, 3, 3], project=proj)
    #     air = op.phase.Air(network=pn)
    #     pickle.dump([pn, air], open('test.pnm', 'wb'))
    #     self.ws.clear()
    #     self.ws.load_project('test.pnm')
    #     self.ws.clear()
    #     os.remove('test.pnm')

    # def test_save_and_load_project_from_pickled_object(self):
    #     a = np.ones((10, ))
    #     pickle.dump(a, open('single_object.pnm', 'wb'))
    #     self.ws.clear()
    #     with pytest.raises(Exception):
    #         self.ws.load_project('single_object.pnm')
    #     b = {'test': a}
    #     pickle.dump(b, open('single_object.pnm', 'wb'))
    #     self.ws.clear()
    #     with pytest.raises(Exception):
    #         self.ws.load_project('single_object.pnm')
    #     os.remove('single_object.pnm')

    # def test_save_and_load_workspace(self):
    #     self.ws.clear()
    #     proj1 = self.ws.new_project('test_proj_1')
    #     proj2 = self.ws.new_project('test_proj_2')
    #     op.network.Cubic(shape=[3, 3, 3], project=proj1, name='net1')
    #     op.network.Cubic(shape=[3, 3, 3], project=proj2, name='net2')
    #     self.ws.save_workspace(filename='workspace_test')
    #     self.ws.clear()
    #     assert 'test_proj_1' not in self.ws.keys()
    #     assert 'test_proj_2' not in self.ws.keys()
    #     self.ws.load_workspace('workspace_test', overwrite=True)
    #     assert 'test_proj_1' in self.ws.keys()
    #     assert 'test_proj_2' in self.ws.keys()
    #     self.ws.clear()
    #     os.remove('workspace_test.pnm')


if __name__ == '__main__':

    t = WorkspaceTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f"Running test {item}")
            t.__getattribute__(item)()
