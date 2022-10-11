import os
import pytest
import numpy as np
import openpnm as op
from pathlib import Path


class ProjectTest:

    def setup_class(self):
        self.ws = op.Workspace()
        self.ws.clear()
        self.net = op.network.Cubic(shape=[2, 2, 2])
        self.proj = self.net.project
        self.phase1 = op.phase.Phase(network=self.net)
        self.phase2 = op.phase.Phase(network=self.net)

    def test_change_project_name_by_assignment(self):
        proj = self.ws.new_project()
        new_name = self.ws._validate_name()
        proj.name = new_name
        assert proj.name == new_name
        assert proj.name in self.ws.keys()

    def test_object_naming(self):
        pn = op.network.Cubic([3, 3, 3], name='bob')
        proj = pn.project
        assert 'bob' in proj.names

    def test_get_locations(self):
        pn = op.network.Cubic([3, 3, 3], name='bob')
        air = op.phase.Air(network=pn)
        assert air.project.get_locations('pore.left').sum() == 9
        assert pn.project.get_locations('pore.left').sum() == 9
        with pytest.raises(KeyError):
            air.project.get_locations('pore.foo')

    # def test_change_simulation_name_by_moving_in_dict(self):
    #     proj = self.ws.new_project()
    #     old_name = proj.name
    #     new_name = self.ws._validate_name()
    #     self.ws[new_name] = proj
    #     assert proj.name == new_name
    #     assert proj.name in self.ws.keys()
    #     assert old_name not in self.ws.keys()

    # def test_grid_access(self):
    #     g = self.proj.grid
    #     with pytest.raises(ValueError):
    #         g.col(self.geo1.name)
    #     r = g.row(self.geo1.name)
    #     assert r == ['geo_01', 'phys_01', 'phys_03']
    #     with pytest.raises(ValueError):
    #         g.row(self.phase1.name)
    #     c = g.col(self.phase1.name)
    #     assert c == ['phase_01', 'phys_01', 'phys_02']
    #     c = g.col(self.phase2.name)
    #     assert c == ['phase_02', 'phys_03', 'phys_04']

    # def test_projectgrid_access(self):
    #     g = self.proj.grid
    #     geoms = g.geometries()
    #     assert geoms == ['geo_01', 'geo_02']
    #     phases = g.phases()
    #     assert phases == ['phase_01', 'phase_02']

    # def test_purge_geom_shallow(self):
    #     proj = self.ws.copy_project(self.net.project)
    #     net = proj.network
    #     geo1 = proj.geometries()['geo_01']
    #     geo2 = proj.geometries()['geo_02']
    #     # Show that geo1 is present to start with
    #     assert 'pore.' + geo1.name in net.keys()
    #     assert 'throat.' + geo1.name in net.keys()
    #     proj.purge_object(geo1)
    #     assert geo1 not in proj
    #     # Ensure geo1 is gone
    #     assert 'pore.' + geo1.name not in net.keys()
    #     assert 'throat.' + geo1.name not in net.keys()
    #     # Check that geo2 is stil there
    #     assert 'pore.' + geo2.name in net.keys()
    #     assert 'throat.' + geo2.name in net.keys()
    #     self.ws.close_project(proj)

    # def test_purge_geom_deep(self):
    #     proj = self.ws.copy_project(self.net.project)
    #     geo1 = proj.geometries()['geo_01']
    #     geo2 = proj.geometries()['geo_02']
    #     phys11 = proj.physics()['phys_01']
    #     phys12 = proj.physics()['phys_02']
    #     phys21 = proj.physics()['phys_03']
    #     phys22 = proj.physics()['phys_04']
    #     proj.purge_object(geo1, deep=True)
    #     assert geo1 not in proj
    #     assert geo2 in proj
    #     assert phys11 not in proj
    #     assert phys12 in proj
    #     assert phys21 not in proj
    #     assert phys22 in proj
    #     self.ws.close_project(proj)

    # def test_purge_phys_shallow(self):
    #     proj = self.ws.copy_project(self.net.project)
    #     phase = proj.phases()['phase_01']
    #     phys1 = proj.physics()['phys_01']
    #     # Show that geo1 is present to start with
    #     assert 'pore.' + phys1.name in phase.keys()
    #     assert 'throat.' + phys1.name in phase.keys()
    #     proj.purge_object(phys1)
    #     assert phys1 not in proj
    #     # Ensure geo1 is gone
    #     assert 'pore.' + phys1.name not in phase.keys()
    #     assert 'throat.' + phys1.name not in phase.keys()
    #     self.ws.close_project(proj)
    #
    # def test_purge_phase_shallow(self):
    #     proj = self.ws.copy_project(self.net.project)
    #     phase = proj.phases()['phase_01']
    #     phys1 = proj.physics()['phys_01']
    #     phys2 = proj.physics()['phys_02']
    #     proj.purge_object(phase)
    #     assert phase not in proj
    #     assert phys1 in proj
    #     assert phys2 in proj
    #     self.ws.close_project(proj)

    # def test_purge_phase_deep(self):
    #     proj = self.ws.copy_project(self.net.project)
    #     phase1 = proj.phases()['phase_01']
    #     phase2 = proj.phases()['phase_02']
    #     phys11 = proj.physics()['phys_01']
    #     phys12 = proj.physics()['phys_02']
    #     phys21 = proj.physics()['phys_03']
    #     phys22 = proj.physics()['phys_04']
    #     proj.purge_object(phase1, deep=True)
    #     assert phase1 not in proj
    #     assert phase2 in proj
    #     assert phys11 not in proj
    #     assert phys12 not in proj
    #     assert phys21 in proj
    #     assert phys22 in proj
    #     self.ws.close_project(proj)

    # def test_purge_network(self):
    #     proj = self.ws.copy_project(self.net.project)
    #     net = proj.network
    #     with pytest.raises(Exception):
    #         proj.purge_object(net)
    #     self.ws.close_project(proj)

    # def test_purge_list(self):
    #     proj = self.ws.copy_project(self.net.project)
    #     phase1 = proj.phases()['phase_01']
    #     phase2 = proj.phases()['phase_02']
    #     phys11 = proj.physics()['phys_01']
    #     phys12 = proj.physics()['phys_02']
    #     phys21 = proj.physics()['phys_03']
    #     phys22 = proj.physics()['phys_04']
    #     proj.purge_object([phase1, phys11, phys12], deep=False)
    #     assert phase1 not in proj
    #     assert phase2 in proj
    #     assert phys11 not in proj
    #     assert phys12 not in proj
    #     assert phys21 in proj
    #     assert phys22 in proj
    #     self.ws.close_project(proj)

    # def test_append_second_network(self):
    #     proj = self.ws.copy_project(self.net.project)
    #     net = proj.network
    #     with pytest.raises(Exception):
    #         proj.append(net)

    # def test_apend_non_openpnm_object(self):
    #     proj = self.ws.copy_project(self.net.project)
    #     with pytest.raises(Exception):
    #         proj.append(1)

    def test_phases(self):
        proj = self.ws.copy_project(self.net.project)
        phases = proj.phases
        # assert 'phase_01' in [p.name for p in phases]
        # assert 'phase_02' in [p.name for p in phases]
        assert len(phases) == 2

    # def test_find_phase_from_physics(self):
    #     proj = self.ws.copy_project(self.net.project)
    #     phys1 = proj.physics()['phys_01']
    #     phys2 = proj.physics()['phys_02']
    #     phase = proj.find_phase(phys1)
    #     assert 'pore.' + phys1.name in phase.keys()
    #     assert 'throat.' + phys1.name in phase.keys()
    #     phase = proj.find_phase(phys1)
    #     assert 'pore.' + phys2.name in phase.keys()
    #     assert 'throat.' + phys2.name in phase.keys()

    # def test_find_phase_from_phase(self):
    #     phases = self.proj.phases
    #     a = self.proj.find_phase(phases[0])
    #     assert a is phases[0]

    # def test_geometries(self):
    #     proj = self.ws.copy_project(self.net.project)
    #     geoms = proj.geometries()
    #     assert 'geo_01' in geoms.keys()
    #     assert 'geo_02' in geoms.keys()
    #     assert len(geoms.keys()) == 2

    # def test_find_geometry(self):
    #     proj = self.ws.copy_project(self.net.project)
    #     phys1 = proj.physics()['phys_01']
    #     phys3 = proj.physics()['phys_03']
    #     geo1 = proj.find_geometry(phys1)
    #     assert geo1.Np == phys1.Np
    #     assert geo1.Nt == phys1.Nt
    #     # Make sure it finds same geometry
    #     assert proj.find_geometry(phys3) is geo1

    # def test_physics(self):
    #     proj = self.ws.copy_project(self.net.project)
    #     physics = proj.physics()
    #     assert 'phys_01' in physics.keys()
    #     assert 'phys_02' in physics.keys()
    #     assert 'phys_03' in physics.keys()
    #     assert 'phys_04' in physics.keys()
    #     assert len(physics.keys()) == 4

    # def test_find_physics_from_geometry(self):
    #     proj = self.proj
    #     geo1 = proj.geometries()['geo_01']
    #     phys1_2 = proj.find_physics(geometry=geo1)
    #     # Ensure two physics were returned
    #     assert len(phys1_2) == 2
    #     # Ensure None is not in list
    #     with pytest.raises(ValueError):
    #         phys1_2.remove(None)
    #     geo2 = proj.geometries()['geo_02']
    #     phys3_4 = proj.find_physics(geometry=geo2)
    #     assert phys3_4 != phys1_2

    # def test_find_physics_from_phase_and_geometry(self):
    #     proj = self.proj
    #     geo1 = proj.geometries()['geo_01']
    #     geo2 = proj.geometries()['geo_02']
    #     phase1 = proj.phases()['phase_01']
    #     phase2 = proj.phases()['phase_02']
    #     phys1 = proj.find_physics(geometry=geo1, phase=phase1)
    #     phys2 = proj.find_physics(geometry=geo2, phase=phase1)
    #     phys3 = proj.find_physics(geometry=geo1, phase=phase2)
    #     phys4 = proj.find_physics(geometry=geo2, phase=phase2)
    #     # Ensure 4 unique physics were returned
    #     assert phys1 is not phys2
    #     assert phys2 is not phys3
    #     assert phys3 is not phys4

    # def test_find_physics_from_phase(self):
    #     proj = self.proj
    #     phase1 = proj.phases()['phase_01']
    #     phase2 = proj.phases()['phase_02']
    #     physics1 = proj.find_physics(phase=phase1)
    #     physics2 = proj.find_physics(phase=phase2)
    #     assert len(physics1) == 2
    #     assert len(physics2) == 2
    #     # Make sure lists are mutually exclusive
    #     assert ~np.all([item in physics2 for item in physics1])

    # def test_find_physics_no_phase_or_geometry(self):
    #     proj = self.proj
    #     a = proj.find_physics()
    #     b = proj.physics().values()
    #     assert np.all([item in b for item in a])

    # def test_find_full_domain_geometry(self):
    #     proj = self.proj
    #     geo1 = proj.geometries()['geo_01']
    #     assert proj.find_full_domain(geo1)._isa() == 'network'

    # def test_find_full_domain_physics(self):
    #     proj = self.proj
    #     phys1 = proj.physics()['phys_01']
    #     assert proj.find_full_domain(phys1)._isa() == 'phase'

    # def test_find_full_domain_phase(self):
    #     proj = self.proj
    #     phase1 = proj.phases()['phase_01']
    #     assert proj.find_full_domain(phase1)._isa() == 'phase'

    # def test_clear(self):
    #     proj = self.ws.copy_project(self.net.project)
    #     assert len(proj) == 3
    #     proj.clear(objtype=['phase'])
    #     assert len(proj) == 7
    #     proj.clear()
    #     assert len(proj) == 0

    # def test_pop(self):
    #     proj = self.ws.copy_project(self.net.project)
    #     geo1 = proj.geometries()['geo_01']
    #     geo2 = proj.pop(1)
    #     assert geo1 is geo2
    #     assert geo1 not in proj

    # def test_insert(self):
    #     proj = self.ws.copy_project(self.net.project)
    #     geo1 = proj.geometries()['geo_01']
    #     geo2 = proj.pop(1)
    #     assert geo1 is geo2
    #     assert geo1 not in proj
    #     proj.insert(1, geo2)

    def test_copy(self):
        proj = self.ws.copy_project(self.net.project)
        proj.network.name = 'foo22'
        proj2 = proj.copy()
        assert proj.name != proj2.name
        assert proj2.network.name == 'foo22'
        assert proj is not proj2
        assert len(proj) == len(proj2)
        assert proj.names == proj2.names

    def test_getitem(self):
        a = self.proj[0]
        b = self.proj[a.name]
        assert a is b

    # def test_save_and_load_object(self):
    #     proj = self.proj
    #     name = proj.network.name
    #     proj.save_object(proj.network)
    #     new_proj = self.ws.new_project()
    #     new_proj.load_object(name+'.net')
    #     assert new_proj.network.name == name
    #     try:
    #         os.remove(name+'.net')
    #     except PermissionError:
    #         print('Could not delete ' + name + '.net')

    # def test_load_object_from_fixture(self):
    #     path = Path(os.path.realpath(__file__),
    #                 '../../../fixtures/OpenPNM-Objects')
    #     filename = Path(path.resolve(), 'net_01.net')
    #     new_proj = self.ws.new_project()
    #     new_proj.load_object(filename)
    #     assert len(new_proj) == 1
    #     assert new_proj.network._isa('network')

    # def test_export_data(self):
    #     fname = 'export_data_tests'
    #     self.proj.export_data(phases=self.phase1, filename=fname,
    #                           filetype='vtk')
    #     os.remove(fname+'.vtp')
    #     self.proj.export_data(phases=self.phase1, filename=fname+'.vtk')
    #     os.remove(fname+'.vtp')
    #     self.proj.export_data(phases=self.phase1, filename=fname,
    #                           filetype='csv')
    #     os.remove(fname+'.csv')
    #     self.proj.export_data(phases=self.phase1, filename=fname,
    #                           filetype='xdmf')
    #     os.remove(fname+'.xmf')
    #     os.remove(fname+'.hdf')
    #     self.proj.export_data(phases=self.phase1, filename=fname,
    #                           filetype='hdf5')
    #     os.remove(fname+'.hdf')
    #     self.proj.export_data(phases=self.phase1, filename=fname,
    #                           filetype='mat')
    #     os.remove(fname+'.mat')

    #     with pytest.raises(Exception):
    #         self.proj.export_data(phases=self.phase1, filename=fname,
    #                               filetype='blah')

    # def test_inspect_pores_and_throats(self):
    #     df = self.proj.inspect_locations(element='pores', indices=[0, 2, 3])
    #     assert df.shape[1] == 3
    #     assert self.net.name + '.pore.coords_X' in df.index
    #     df = self.proj.inspect_locations(element='throats', indices=[0, 1, 3])
    #     assert df.shape[1] == 3
    #     assert self.net.name + '.throat.conns_head' in df.index


if __name__ == '__main__':

    t = ProjectTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            t.__getattribute__(item)()
