import OpenPNM
import numpy as np


class GenericAlgorithmTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        self.phase = OpenPNM.Phases.GenericPhase(network=self.net)
        self.alg = OpenPNM.Algorithms.GenericLinearTransport(network=self.net,
                                                             phase=self.phase)
