import pytest
import OpenPNM


class GenericPhaseTest:
    def setup_class(self):
        self.net = OpenPNM.Network.Cubic(shape=[5, 5, 5])

    def test_init_w_no_network(self):
        OpenPNM.Phases.GenericPhase()

    def test_init_w_components(self):
        comp1 = OpenPNM.Phases.GenericPhase(network=self.net)
        comp2 = OpenPNM.Phases.GenericPhase(network=self.net)
        OpenPNM.Phases.GenericPhase(network=self.net,
                                    components=[comp1, comp2])

    def test_set_component_add(self):
        comp1 = OpenPNM.Phases.GenericPhase(network=self.net)
        comp2 = OpenPNM.Phases.GenericPhase(network=self.net)
        phase = OpenPNM.Phases.GenericPhase(network=self.net)
        phase.set_component(comp1)
        phase.set_component(comp2)

    def test_set_component_add_twice(self):
        comp1 = OpenPNM.Phases.GenericPhase(network=self.net)
        phase = OpenPNM.Phases.GenericPhase(network=self.net)
        phase.set_component(comp1)
        with pytest.raises(Exception):
            phase.set_components(comp1)

    def test_set_component_remove(self):
        comp1 = OpenPNM.Phases.GenericPhase(network=self.net)
        comp2 = OpenPNM.Phases.GenericPhase(network=self.net)
        phase = OpenPNM.Phases.GenericPhase(network=self.net,
                                            components=[comp1, comp2])
        phase.set_component(comp1, mode='remove')
        phase.set_component(comp2, mode='remove')

    def test_set_component_remove_twice(self):
        comp1 = OpenPNM.Phases.GenericPhase(network=self.net)
        comp2 = OpenPNM.Phases.GenericPhase(network=self.net)
        phase = OpenPNM.Phases.GenericPhase(network=self.net,
                                            components=[comp1, comp2])
        phase.set_component(comp1, mode='remove')
        with pytest.raises(Exception):
            phase.set_component(comp1, mode='remove')
