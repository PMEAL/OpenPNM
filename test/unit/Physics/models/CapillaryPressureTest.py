import OpenPNM


class CapillaryPressureTest:

    def setup_class(self):
        self.pn = OpenPNM.Network.TestNet()
        self.phase = OpenPNM.Phases.TestPhase()
        self.physics = OpenPNM.Physics.TestPhysics()

    def test_washburn(self):
        pass
