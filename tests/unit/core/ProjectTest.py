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

    def test_purge_geom(self):
        proj = self.ws.copy_project(self.net.project)
        net = proj.network
        geo1 = proj.geometries()['geo_01']
        geo2 = proj.geometries()['geo_02']
        # Show that geo1 is present to start with
        assert 'pore.' + geo1.name in net.keys()
        assert 'throat.' + geo1.name in net.keys()
        proj.purge_object(geo1)
        assert geo1 not in proj
        # Ensure geo1 is gone
        assert 'pore.' + geo1.name not in net.keys()
        assert 'throat.' + geo1.name not in net.keys()
        # Check that geo2 is stil there
        assert 'pore.' + geo2.name in net.keys()
        assert 'throat.' + geo2.name in net.keys()
        self.ws.close_project(proj)

    def test_purge_phys(self):
        proj = self.ws.copy_project(self.net.project)
        phase = proj.phases()['phase_01']
        phys1 = proj.physics()['phys_01']
        # Show that geo1 is present to start with
        assert 'pore.' + phys1.name in phase.keys()
        assert 'throat.' + phys1.name in phase.keys()
        proj.purge_object(phys1)
        assert phys1 not in proj
        # Ensure geo1 is gone
        assert 'pore.' + phys1.name not in phase.keys()
        assert 'throat.' + phys1.name not in phase.keys()
        self.ws.close_project(proj)

    def test_purge_phase(self):
        proj = self.ws.copy_project(self.net.project)
        phase = proj.phases()['phase_01']
        phys1 = proj.physics()['phys_01']
        phys2 = proj.physics()['phys_02']
        proj.purge_object(phase)
        assert phys1 not in proj
        assert phys2 not in proj
        assert phase not in proj
        self.ws.close_project(proj)

    def test_purge_network(self):
        proj = self.ws.copy_project(self.net.project)
        net = proj.network
        with pytest.raises(Exception):
            proj.purge_object(net)
        self.ws.close_project(proj)

    def test_append_second_network(self):
        proj = self.ws.copy_project(self.net.project)
        net = proj.network
        with pytest.raises(Exception):
            proj.append(net)

    def test_apend_non_openpnm_object(self):
        proj = self.ws.copy_project(self.net.project)
        with pytest.raises(Exception):
            proj.append(1)

    def test_phases(self):
        proj = self.ws.copy_project(self.net.project)
        phases = proj.phases()
        assert 'phase_01' in phases.keys()
        assert 'phase_02' in phases.keys()
        assert len(phases.keys()) == 2

    def test_find_phase(self):
        proj = self.ws.copy_project(self.net.project)
        phys1 = proj.physics()['phys_01']
        phys2 = proj.physics()['phys_02']
        phase = proj.find_phase(phys1)
        assert 'pore.' + phys1.name in phase.keys()
        assert 'throat.' + phys1.name in phase.keys()
        phase = proj.find_phase(phys1)
        assert 'pore.' + phys2.name in phase.keys()
        assert 'throat.' + phys2.name in phase.keys()

    def test_geometries(self):
        proj = self.ws.copy_project(self.net.project)
        geoms = proj.geometries()
        assert 'geo_01' in geoms.keys()
        assert 'geo_02' in geoms.keys()
        assert len(geoms.keys()) == 2

    def test_find_geometry(self):
        proj = self.ws.copy_project(self.net.project)
        phys1 = proj.physics()['phys_01']
        phys3 = proj.physics()['phys_03']
        geo1 = proj.find_geometry(phys1)
        assert geo1.Np == phys1.Np
        assert geo1.Nt == phys1.Nt
        # Make sure it finds same geometry
        assert proj.find_geometry(phys3) is geo1

    def test_physics(self):
        proj = self.ws.copy_project(self.net.project)
        physics = proj.physics()
        assert 'phys_01' in physics.keys()
        assert 'phys_02' in physics.keys()
        assert 'phys_03' in physics.keys()
        assert 'phys_04' in physics.keys()
        assert len(physics.keys()) == 4

    def test_find_physics(self):
        proj = self.ws.copy_project(self.net.project)
        geo1 = proj.geometries()['geo_01']
        phys1_2 = proj.find_physics(geometry=geo1)
        # Ensure two physics were returned
        assert len(phys1_2) == 2
        # Ensure None is not in list
        with pytest.raises(ValueError):
            phys1_2.remove(None)
        geo2 = proj.geometries()['geo_02']
        phys3_4 = proj.find_physics(geometry=geo2)
        assert phys3_4 != phys1_2
        phase1 = proj.phases()['phase_01']
        phase2 = proj.phases()['phase_02']
        phys1 = proj.find_physics(geometry=geo1, phase=phase1)
        phys2 = proj.find_physics(geometry=geo2, phase=phase1)
        phys3 = proj.find_physics(geometry=geo1, phase=phase2)
        phys4 = proj.find_physics(geometry=geo2, phase=phase2)
        # Ensure 4 unique physics were returned
        assert phys1 is not phys2
        assert phys2 is not phys3
        assert phys3 is not phys4

    def test_clear(self):
        proj = self.ws.copy_project(self.net.project)
        assert len(proj) == 9
        proj.clear(objtype=['phase'])
        assert len(proj) == 3
        proj.clear()
        assert len(proj) == 0



if __name__ == '__main__':

    t = ProjectTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
