import OpenPNM


class ElectricalConductivityTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[3, 3, 3])
        self.phase = OpenPNM.Phases.GenericPhase(network=self.net)
        self.phase['pore.intrinsic_conductivity'] = 1
        self.phase['pore.volume_fraction'] = 0.5

    def test_percolating_continua(self):
        f = OpenPNM.Phases.models.electrical_conductivity.percolating_continua
        self.phase.models.add(propname='pore.effective_conductivity',
                              model=f,
                              phi_crit=0.25,
                              tau=2)
