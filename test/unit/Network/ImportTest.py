import OpenPNM
import scipy as sp
from os.path import join


class Import:
    def setup_class(self):
        fname = join(FIXTURE_DIR, 'import_data.csv')
        self.net = OpenPNM.Network.Import(file=fname)
