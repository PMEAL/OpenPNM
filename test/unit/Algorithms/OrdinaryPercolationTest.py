import OpenPNM as op
import scipy as sp
mgr = op.Base.Workspace()
mgr.loglevel = 60


class OrdinaryPercolationTest:
    def setup_test(self):
        self.net = op.Network.Cubic(shape=[5, 5, 5])
        self.geo = op.Geometry.Toray090(network=self.net,
                                        pores=self.net.Ps,
                                        throats=self.net.Ts)
        self.phase = op.Phases.Water(network=self.net)
        self.phys = op.Physics.Standard(network=self.net,
                                        pores=self.net.Ps,
                                        throats=self.net.Ts)
        self.OP = op.Algorithms.OrdinaryPercolation(network=self.net,
                                                    invading_phase=self.phase)
        Ps = self.net.pores(labels=['bottom_boundary'])
        self.OP.run(inlets=Ps)
        self.OP.return_results(Pc=7000)
        lpf = self.OP.evaluate_late_pore_filling(Pc=8000)
        assert sp.size(lpf) == self.net.Np
