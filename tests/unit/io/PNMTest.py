import py
import os
import shutil
import pytest
import numpy as np
import openpnm as op
ws = op.Workspace()


class PNMTest:

    def setup_class(self):
        ws.settings['local_data'] = True

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_save_and_reload(self):
        f = 'test1.pnm'
        pn = op.network.Cubic(shape=[3, 3, 3])
        pn['pore.random'] = np.random.rand(pn.Np)
        op.io.PNM.save_project(project=pn.project, filename=f)
        ws.clear()
        proj = op.io.PNM.load_project(f)
        net = proj.network
        assert np.all(net['pore.random'] == pn['pore.random'])
        shutil.rmtree(f, ignore_errors=True)
        os.remove("test1.pnm")

    def test_save_and_load_with_models(self):
        f = 'test2.pnm'
        pn = op.network.Cubic(shape=[3, 3, 3])
        op.io.PNM.save_project(project=pn.project, filename=f)
        ws.clear()
        proj = op.io.PNM.load_project(f)
        net = proj.network
        net.regenerate_models('pore.coordination_number')
        assert 'pore.coordination_number' not in pn.keys()
        assert 'pore.coordination_number' in net.keys()
        shutil.rmtree(f, ignore_errors=True)
        os.remove("test2.pnm")

    def test_save_and_load_with_local_custom_model(self):
        f = 'test3.pnm'

        def test(target):
            return 1.0

        pn = op.network.Cubic(shape=[3, 3, 3])
        pn.add_model(propname='pore.test', model=test)
        assert np.all(pn['pore.test'] == 1.0)
        op.io.PNM.save_project(project=pn.project, filename=f)
        ws.clear()
        proj = op.io.PNM.load_project(f)
        net = proj.network
        assert 'pore.test' in net.models.keys()
        assert np.all(net['pore.test'] == 1.0)
        # Now delete data and re-run
        del net['pore.test']
        net.regenerate_models()
        with pytest.raises(Exception):
            assert np.all(net['pore.test'] == 1.0)
        shutil.rmtree(f, ignore_errors=True)
        os.remove("test3.pnm")

    def test_save_and_load_with_imported_custom_model(self):
        from custom_code import test
        f = 'test4.pnm'
        pn = op.network.Cubic(shape=[3, 3, 3])
        pn.add_model(propname='pore.test', model=test)
        assert np.all(pn['pore.test'] == 1.0)
        op.io.PNM.save_project(project=pn.project, filename=f)
        ws.clear()
        proj = op.io.PNM.load_project(f)
        net = proj.network
        assert 'pore.test' in net.models.keys()
        assert np.all(net['pore.test'] == 1.0)
        # Now delete data and re-run
        del net['pore.test']
        net.regenerate_models()
        assert np.all(net['pore.test'] == 1.0)
        shutil.rmtree(f, ignore_errors=True)
        os.remove("test4.pnm")


if __name__ == '__main__':
    # All the tests in this file can be run with 'playing' this file
    t = PNMTest()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=py.path.local())
