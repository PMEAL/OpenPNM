import os
import openpnm as op
from pathlib import Path
import pickle


class PoreSpyTest:

    def setup_class(self):
        self.path = os.path.join(Path(__file__).parent.resolve(), "berea.net")
        with open(self.path, 'rb') as f:
            self.net = pickle.load(f)

    def test_load_PoreSpy_from_pickle(self):
        proj = op.io.network_from_porespy(self.net)
        net = proj.network
        assert net.Np == 1637
        assert net.Nt == 2785

    def test_load_PoreSpy_from_file(self):
        proj = op.io.network_from_porespy(filename=self.path)
        net = proj.network
        assert net.Np == 1637
        assert net.Nt == 2785

    def test_param_keys_are_routed_to_params(self):
        net_with_params = dict(self.net)
        net_with_params['param.voxel_size'] = 1.5e-6
        net_with_params['param.ndim'] = 3
        proj = op.io.network_from_porespy(net_with_params)
        net = proj.network
        assert net['param.voxel_size'] == 1.5e-6
        assert net['param.ndim'] == 3
        # Ensure they survive a project save/load round-trip
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            fname = os.path.join(tmp, 'proj')
            op.Workspace().save_project(project=net.project, filename=fname)
            ws = op.Workspace()
            ws.clear()
            ws.load_project(filename=fname + '.pnm')
            loaded = list(ws.values())[0][0]
            assert loaded['param.voxel_size'] == 1.5e-6
            assert loaded['param.ndim'] == 3


if __name__ == '__main__':
    # All the tests in this file can be run with 'playing' this file
    t = PoreSpyTest()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
