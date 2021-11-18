import os
import py
import pytest
import pickle
import openpnm as op


class PickleTest:

    def setup_class(self):
        pass

    def teardown_class(self, tmpdir=None):
        os.remove('test.pkl')
        os.remove("proj.pkl")
        os.remove("pn.pkl")

    def test_create_save_and_load_workspace(self, tmpdir=None):
        ws = op.Workspace()
        ws.clear()
        pn = op.network.Cubic(shape=[3, 3, 3])
        assert len(ws) == 1
        op.io.Pickle.save_workspace('test.pkl')
        ws = op.io.Pickle.load_workspace('test.pkl', overwrite=False)
        assert isinstance(ws, dict)
        assert len(ws) == 2
        ws = op.io.Pickle.load_workspace('test.pkl', overwrite=True)
        assert len(ws) == 1

    def test_create_save_and_load_project(self, tmpdir=None):
        ws = op.Workspace()
        ws.clear()
        pn = op.network.Cubic(shape=[3, 3, 3])
        assert len(ws) == 1
        op.io.Pickle.save_project(project=pn.project, filename='proj.pkl')
        proj = op.io.Pickle.load_project('proj.pkl')
        assert isinstance(proj, list)
        assert len(ws) == 2

    def test_load_project_with_multiple_projects(self, tmpdir=None):
        ws = op.Workspace()
        ws.clear()
        pn1 = op.network.Cubic(shape=[3, 3, 3])
        pn2 = op.network.Cubic(shape=[3, 3, 3])
        assert len(ws) == 2
        op.io.Pickle.save_workspace(filename='proj.pkl')
        with pytest.raises(Exception):
            proj = op.io.Pickle.load_project('proj.pkl')

    def test_load_handmade_project(self, tmpdir=None):
        ws = op.Workspace()
        ws.clear()
        pn = op.network.Cubic(shape=[3, 3, 3])
        new_proj = [pn]
        pickle.dump(new_proj, open('proj.pkl', 'wb'))
        proj = op.io.Pickle.load_project('proj.pkl')
        assert isinstance(proj, op.Project)
        assert len(ws.keys()) == 2

    def test_load_poorly_made_project(self, tmpdir=None):
        proj = 5
        pickle.dump(proj, open('proj.pkl', 'wb'))
        with pytest.raises(Exception):
            proj = op.io.Pickle.load_project('proj.pkl')

    def test_load_from_saved_dict(self, tmpdir=None):
        ws = op.Workspace()
        ws.clear()
        pn = op.network.Cubic(shape=[3, 3, 3])
        op.io.Pickle.save_workspace(filename='test.pkl')
        ws = op.io.Pickle.load_workspace(filename='test.pkl',
                                            overwrite=True)
        assert len(ws.keys()) == 1
        ws = op.io.Pickle.load_workspace(filename='test.pkl',
                                            overwrite=False)
        assert len(ws.keys()) == 2
        assert isinstance(ws, op.Workspace)

    def test_load_workspace_from_poorly_made_dict(self):
        ws = op.Workspace()
        ws.clear()
        pn = op.network.Cubic(shape=[3, 3, 3])
        pickle.dump(pn, open('pn.pkl', 'wb'))
        with pytest.raises(Exception):
            ws = op.io.Pickle.load_workspace('pn.pkl')


if __name__ == '__main__':
    # All the tests in this file can be run with 'playing' this file
    t = PickleTest()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=py.path.local())
    t.teardown_class()
