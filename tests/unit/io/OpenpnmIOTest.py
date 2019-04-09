import openpnm as op
import pytest
import py
import os
import pickle


class OpenpnmIOTest:

    def setup_class(self):
        pass

    def test_create_save_and_load_workspace(self, tmpdir=None):
        ws = op.Workspace()
        ws.clear()
        pn = op.network.Cubic(shape=[3, 3, 3])
        assert len(ws) == 1
        op.io.OpenpnmIO.save_workspace('test.pnm')
        ws = op.io.OpenpnmIO.load_workspace('test.pnm', overwrite=False)
        assert isinstance(ws, dict)
        assert len(ws) == 2
        ws = op.io.OpenpnmIO.load_workspace('test.pnm', overwrite=True)
        assert len(ws) == 1
        os.remove('test.pnm')

    def test_create_save_and_load_project(self, tmpdir=None):
        ws = op.Workspace()
        ws.clear()
        pn = op.network.Cubic(shape=[3, 3, 3])
        assert len(ws) == 1
        op.io.OpenpnmIO.save_project(project=pn.project, filename='proj.pnm')
        proj = op.io.OpenpnmIO.load_project('proj.pnm')
        assert isinstance(proj, list)
        assert len(ws) == 2
        os.remove('proj.pnm')

    def test_load_project_with_multiple_projects(self, tmpdir=None):
        ws = op.Workspace()
        ws.clear()
        pn1 = op.network.Cubic(shape=[3, 3, 3])
        pn2 = op.network.Cubic(shape=[3, 3, 3])
        assert len(ws) == 2
        op.io.OpenpnmIO.save_workspace(filename='proj.pnm')
        with pytest.raises(Exception):
            proj = op.io.OpenpnmIO.load_project('proj.pnm')
        os.remove('proj.pnm')

    def test_load_handmade_project(self, tmpdir=None):
        ws = op.Workspace()
        ws.clear()
        pn = op.network.Cubic(shape=[3, 3, 3])
        new_proj = [pn]
        pickle.dump(new_proj, open('proj.pnm', 'wb'))
        proj = op.io.OpenpnmIO.load_project('proj.pnm')
        assert isinstance(proj, op.Project)
        assert len(ws.keys()) == 2
        os.remove('proj.pnm')

    def test_load_poorly_made_project(self, tmpdir=None):
        proj = 5
        pickle.dump(proj, open('proj.pnm', 'wb'))
        with pytest.raises(Exception):
            proj = op.io.OpenpnmIO.load_project('proj.pnm')
        os.remove('proj.pnm')

    def test_load_from_saved_dict(self, tmpdir=None):
        ws = op.Workspace()
        ws.clear()
        pn = op.network.Cubic(shape=[3, 3, 3])
        op.io.OpenpnmIO.save_workspace(filename='test.pnm')
        ws = op.io.OpenpnmIO.load(filename='test.pnm')
        assert len(ws.keys()) == 1
        assert isinstance(ws, op.Workspace)

    def test_load_workspace_from_poorly_made_dict(self):
        ws = op.Workspace()
        ws.clear()
        pn = op.network.Cubic(shape=[3, 3, 3])
        pickle.dump(pn, open('pn.pnm', 'wb'))
        with pytest.raises(Exception):
            ws = op.io.OpenpnmIO.load_workspace('pn.pnm')


if __name__ == '__main__':
    # All the tests in this file can be run with 'playing' this file
    t = OpenpnmIOTest()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=py.path.local())
