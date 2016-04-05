import OpenPNM
import pytest


class GenericAlgorithmTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        self.phase = OpenPNM.Phases.GenericPhase(network=self.net)
        self.alg = OpenPNM.Algorithms.GenericAlgorithm(network=self.net,
                                                       phase=self.phase)

    def test_non_implemented_methods(self):
        with pytest.raises(NotImplementedError):
            self.alg.run()
        with pytest.raises(NotImplementedError):
            self.alg.set_boundary_conditions()
        with pytest.raises(NotImplementedError):
            self.alg.return_results()
