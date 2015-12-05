import OpenPNM
import scipy as sp


class StandardTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])
        self.geo = OpenPNM.Geometry.Stick_and_Ball(network=self.net,
                                                   pores=self.net.Ps,
                                                   throats=self.net.Ts)
        self.phase = OpenPNM.Phases.Air(network=self.net)
        self.phys = OpenPNM.Physics.Standard(network=self.net,
                                             geometry=self.geo,
                                             phase=self.phase)

    def test_reassign_phase(self):
        phase = OpenPNM.Phases.Water(network=self.net)
        # Store current values from phys for later
        phys_vals = self.phys['throat.hydraulic_conductance']
        # Confirm that phase is not associated with
        flag = False
        try:
            vals = phase['throat.hydraulic_conductance']
        except:
            flag = True
        assert flag
        # Now reassign phys to phase
        self.phys.parent_phase = phase
        # And recheck for throat conductance on phase
        try:
            vals = phase['throat.hydraulic_conductance']
        except:
            flag = True
        assert flag




