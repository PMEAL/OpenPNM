import OpenPNM
import scipy as sp
from os.path import join


class ImportTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Import()
        ctrl = OpenPNM.Base.Controller()
        ctrl.loglevel = 50

    def test_from_csv(self):
        fname = join(FIXTURE_DIR, 'test_load_csv_no_phases.csv')
        net = OpenPNM.Network.Import()
        assert sorted(net.keys()) == ['pore.all', 'throat.all']
        net.from_csv(filename=fname, mode='overwrite')
        assert net.Np == 27
        assert 'pore.seed' in net.keys()
        assert sp.sum(net['pore.seed']) > 0
        net['pore.seed'] = 0
        assert sp.sum(net['pore.seed']) == 0
        net.from_csv(filename=fname, mode='overwrite')
        assert sp.sum(net['pore.seed']) > 0
        net['pore.seed'] = 0
        net.from_csv(filename=fname, mode='add')
        assert sp.sum(net['pore.seed']) == 0
        del net['pore.seed']
        net.from_csv(filename=fname, mode='add')
        assert sp.sum(net['pore.seed']) > 0

    def test_from_mat(self):
        fname = join(FIXTURE_DIR, 'test_load_mat_no_phases.mat')
        net = OpenPNM.Network.Import()
        assert sorted(net.keys()) == ['pore.all', 'throat.all']
        net.from_mat(filename=fname, mode='overwrite')
        assert net.Np == 27
        assert 'pore.index' in net.keys()
        assert sp.sum(net['pore.index']) > 0
        net['pore.index'] = 0
        assert sp.sum(net['pore.index']) == 0
        net.from_mat(filename=fname, mode='overwrite')
        assert sp.sum(net['pore.index']) > 0
        net['pore.index'] = 0
        net.from_mat(filename=fname, mode='add')
        assert sp.sum(net['pore.index']) == 0
        del net['pore.index']
        net.from_mat(filename=fname, mode='add')
        assert sp.sum(net['pore.index']) > 0

    def test_from_vtk(self):
        fname = join(FIXTURE_DIR, 'test_load_vtk_no_phases.vtp')
        net = OpenPNM.Network.Import()
        assert sorted(net.keys()) == ['pore.all', 'throat.all']
        net.from_vtk(filename=fname, mode='overwrite')
        assert net.Np == 27
        assert 'pore.seed' in net.keys()
        assert sp.sum(net['pore.seed']) > 0
        net['pore.seed'] = 0
        assert sp.sum(net['pore.seed']) == 0
        net.from_vtk(filename=fname, mode='overwrite')
        assert sp.sum(net['pore.seed']) > 0
        net['pore.seed'] = 0
        net.from_vtk(filename=fname, mode='add')
        assert sp.sum(net['pore.seed']) == 0
        del net['pore.seed']
        net.from_vtk(filename=fname, mode='add')
        assert sp.sum(net['pore.seed']) > 0

    def test_from_yaml(self):
        fname = join(FIXTURE_DIR, 'test_load_yaml.yaml')
        net = OpenPNM.Network.Import()
        assert sorted(net.keys()) == ['pore.all', 'throat.all']
        net.from_yaml(filename=fname, mode='overwrite')
        assert net.Np == 9
        assert net.Nt == 12
        assert sp.shape(net['pore.coords']) == (9, 3)
        assert sp.shape(net['throat.conns']) == (12, 2)
        assert 'pore.diameter' in net.keys()
