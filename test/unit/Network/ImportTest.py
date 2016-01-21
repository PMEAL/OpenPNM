import OpenPNM
import scipy as sp
from os.path import join


class ImportTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Import()
        FIXTURE_DIR = 'C:\\Users\\Jeff\\Dropbox\\Flash Sync\\Code\\Git\\OpenPNM\\test\\fixtures'

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
