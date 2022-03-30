import openpnm as op


class ElectricalConductivityTest:
    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.phase = op.phase.GenericPhase(network=self.net)
        self.phase['pore.intrinsic_conductivity'] = 1
        self.phase['pore.volume_fraction'] = 0.5

    def test_percolating_continua(self):
        f = op.models.phase.electrical_conductivity.percolating_continua
        self.phase.add_model(propname='pore.effective_conductivity',
                             model=f,
                             phi_crit=0.25,
                             tau=2)


if __name__ == '__main__':

    t = ElectricalConductivityTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
