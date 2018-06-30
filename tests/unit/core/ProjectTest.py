import openpnm as op
import scipy as sp
import pytest
from pathlib import Path
import os


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

    def test_find_phase_from_physics(self):
        proj = self.ws.copy_project(self.net.project)
        phys1 = proj.physics()['phys_01']
        phys2 = proj.physics()['phys_02']
        phase = proj.find_phase(phys1)
        assert 'pore.' + phys1.name in phase.keys()
        assert 'throat.' + phys1.name in phase.keys()
        phase = proj.find_phase(phys1)
        assert 'pore.' + phys2.name in phase.keys()
        assert 'throat.' + phys2.name in phase.keys()

    def test_find_phase_from_phase(self):
        phases = list(self.proj.phases().values())
        a = self.proj.find_phase(phases[0])
        assert a is phases[0]

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

    def test_find_physics_from_geometry(self):
        proj = self.proj
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

    def test_find_physics_from_phase_and_geometry(self):
        proj = self.proj
        geo1 = proj.geometries()['geo_01']
        geo2 = proj.geometries()['geo_02']
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

    def test_find_physics_from_phase(self):
        proj = self.proj
        phase1 = proj.phases()['phase_01']
        phase2 = proj.phases()['phase_02']
        physics1 = proj.find_physics(phase=phase1)
        physics2 = proj.find_physics(phase=phase2)
        assert len(physics1) == 2
        assert len(physics2) == 2
        # Make sure lists are mutually exclusive
        assert ~sp.all([item in physics2 for item in physics1])

    def test_find_physics_no_phase_or_geometry(self):
        proj = self.proj
        a = proj.find_physics()
        b = proj.physics().values()
        assert sp.all([item in b for item in a])

    def test_clear(self):
        proj = self.ws.copy_project(self.net.project)
        assert len(proj) == 9
        proj.clear(objtype=['phase'])
        assert len(proj) == 3
        proj.clear()
        assert len(proj) == 0

    def test_getitem(self):
        a = self.proj[0]
        b = self.proj[a.name]
        assert a is b

    def test_comments(self):
        proj = self.proj
        proj.comments = 'test comment'
        assert 'test comment' in proj._comments.values()
        proj.comments

    def test_print(self):
        proj = self.proj
        s = proj.__str__()
        # 13 rows
        assert len(s.split('\n')) == 13

    def test_save_and_load_object(self):
        proj = self.proj
        name = proj.network.name
        proj.save_object(proj.network)
        new_proj = self.ws.new_project()
        new_proj.load_object(name+'.net')
        assert new_proj.network.name == name
        os.remove(name+'.net')

    def test_load_object_from_fixture(self):
        path = Path(os.path.realpath(__file__),
                    '../../../fixtures/OpenPNM-Objects')
        filename = Path(path.resolve(), 'net_01.net')
        new_proj = self.ws.new_project()
        new_proj.load_object(filename)
        assert len(new_proj) == 1
        assert new_proj.network._isa('network')

    def test_dump_and_fetch_data(self):
        proj = self.ws.copy_project(self.proj)
        proj._dump_data()
        # Ensure only pore.coords and throat.conns are found
        assert sum([len(item.props()) for item in proj]) == 2
        proj._fetch_data()
        assert sp.any([len(item.props()) for item in proj])
        os.remove(proj.name+'.hdf5')

    def test_export_data(self):
        fname = 'export_data_tests'
        self.proj.export_data(network=self.net, phases=self.phase1,
                              filename=fname, filetype='vtp')
        os.remove(fname+'.vtp')
        self.proj.export_data(network=self.net, phases=self.phase1,
                              filename=fname+'.vtp')
        os.remove(fname+'.vtp')
        self.proj.export_data(network=self.net, phases=self.phase1,
                              filename=fname, filetype='csv')
        os.remove(fname+'.csv')
#        self.proj.export_data(network=self.net, phases=self.phase1,
#                              filename=fname, filetype='xmf')
#        os.remove(fname+'.xmf')
#        os.remove(fname+'.hdf')
        self.proj.export_data(network=self.net, phases=self.phase1,
                              filename=fname, filetype='hdf')
        os.remove(fname+'.hdf')
        self.proj.export_data(network=self.net, phases=self.phase1,
                              filename=fname, filetype='mat')
        os.remove(fname+'.mat')


if __name__ == '__main__':

    t = ProjectTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
