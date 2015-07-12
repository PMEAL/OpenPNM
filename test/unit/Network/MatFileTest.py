import OpenPNM
import scipy as sp
from os.path import join


class MatFileTest:
    def setup_class(self):
        fname = join(FIXTURE_DIR, 'example_network.mat')
        self.net = OpenPNM.Network.MatFile(filename=fname)
