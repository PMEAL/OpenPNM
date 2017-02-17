import OpenPNM as op
import scipy as sp
mgr = op.Base.Workspace()
mgr.loglevel = 60


class OrdinaryPercolationTest:

    def setup_class(self):
        self.net = op.Network.Cubic(shape=[5, 5, 5])
        self.geo = op.Geometry.Toray090(network=self.net,
                                        pores=self.net.Ps,
                                        throats=self.net.Ts)
        self.phase = op.Phases.Water(network=self.net)
        self.phys = op.Physics.Standard(network=self.net,
                                        phase=self.phase,
                                        pores=self.net.Ps,
                                        throats=self.net.Ts)
        self.OP1 = op.Algorithms.OrdinaryPercolation(network=self.net,
                                                     invading_phase=self.phase)
        Ps = self.net.pores(labels=['bottom'])
        self.OP1.run(inlets=Ps)
        self.OP1.return_results(Pc=7000)
        lpf = self.OP1.evaluate_late_pore_filling(Pc=8000)
        assert sp.size(lpf) == self.net.Np

    def test_site_percolation(self):
        self.OP2 = op.Algorithms.OrdinaryPercolation(network=self.net,
                                                     invading_phase=self.phase,
                                                     percolation_type='site')
        Ps = self.net.pores(labels=['bottom'])
        self.OP2.run(inlets=Ps)
        self.OP2.return_results(Pc=7000)
        lpf = self.OP2.evaluate_late_pore_filling(Pc=8000)
        assert sp.size(lpf) == self.net.Np
