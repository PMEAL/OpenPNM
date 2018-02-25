import openpnm as op
import scipy as sp
import pytest


class ProjectTest:

    def setup_class(self):
        self.ws = op.core.Workspace()
        self.ws.clear()
        self.proj = self.ws.new_project()
        self.net = op.network.Cubic(shape=[2, 2, 2], project=self.proj)
        Ps = self.net.pores('top')
        Ts = self.net.find_neighbor_throats(pores=Ps)
        self.geo1 = op.geometry.GenericGeometry(network=self.net, pores=Ps,
                                                throats=Ts)
        Ps = self.net.pores('bottom')
        Ts = ~self.net.tomask(throats=Ts)
        self.geo2 = op.geometry.GenericGeometry(network=self.net, pores=Ps,
                                                throats=Ts)
        self.phase1 = op.phases.GenericPhase(network=self.net)
        self.phase2 = op.phases.GenericPhase(network=self.net)
        self.phys11 = op.physics.GenericPhysics(network=self.net,
                                                phase=self.phase1,
                                                geometry=self.geo1)
        self.phys12 = op.physics.GenericPhysics(network=self.net,
                                                phase=self.phase1,
                                                geometry=self.geo2)
        self.phys21 = op.physics.GenericPhysics(network=self.net,
                                                phase=self.phase2,
                                                geometry=self.geo1)
        self.phys22 = op.physics.GenericPhysics(network=self.net,
                                                phase=self.phase2,
                                                geometry=self.geo2)

    def test_change_simulation_name_by_assignment(self):
        proj = self.ws.new_project()
        new_name = self.ws._gen_name()
        proj.name = new_name
        assert proj.name == new_name
        assert proj.name in self.ws.keys()

    def test_change_simulation_name_by_moving_in_dict(self):
        proj = self.ws.new_project()
        old_name = proj.name
        new_name = self.ws._gen_name()
        self.ws[new_name] = proj
        assert proj.name == new_name
        assert proj.name in self.ws.keys()
        assert old_name not in self.ws.keys()

    def test_grid_printing(self):
        d = self.proj.grid
        assert d == {'geo_01': {'phase_01': 'phys_01', 'phase_02': 'phys_03'},
                     'geo_02': {'phase_01': 'phys_02', 'phase_02': 'phys_04'}}

        s = "――――――――――――――――――――――――――――――――――――――――――――――――\n" + \
            "|     net_01   |    phase_01   |    phase_02   |\n" + \
            "――――――――――――――――――――――――――――――――――――――――――――――――\n" + \
            "|     geo_01   |    phys_01    |    phys_03    |\n" + \
            "――――――――――――――――――――――――――――――――――――――――――――――――\n" + \
            "|     geo_02   |    phys_02    |    phys_04    |\n" + \
            "――――――――――――――――――――――――――――――――――――――――――――――――"
        assert print(self.proj.grid) == print(s)


if __name__ == '__main__':

    t = ProjectTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
