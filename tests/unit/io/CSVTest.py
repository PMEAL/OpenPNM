import openpnm as op
import pytest
import py
import os


class CSVTest:

    def setup_class(self):
        ws = op.Workspace()
        ws.settings['local_data'] = True
        self.net = op.network.Cubic(shape=[2, 2, 2])
        Ps = [0, 1, 2, 3]
        Ts = self.net.find_neighbor_throats(pores=Ps)
        self.geo_1 = op.geometry.GenericGeometry(network=self.net,
                                                 pores=Ps, throats=Ts)
        self.geo_1['pore.boo'] = 1
        self.geo_1['throat.boo'] = 1
        Ps = [4, 5, 6, 7]
        Ts = self.net.find_neighbor_throats(pores=Ps, mode='xnor')
        self.geo_2 = op.geometry.GenericGeometry(network=self.net,
                                                 pores=Ps, throats=Ts)
        self.geo_2['pore.boo'] = 1
        self.geo_2['throat.boo'] = 1

        self.phase_1 = op.phases.GenericPhase(network=self.net)
        self.phase_1['pore.bar'] = 2
        self.phase_1['throat.bar'] = 2
        self.phase_2 = op.phases.GenericPhase(network=self.net)
        self.phase_2['pore.bar'] = 2
        self.phase_2['throat.bar'] = 2

        self.phys_1 = op.physics.GenericPhysics(network=self.net,
                                                phase=self.phase_1,
                                                geometry=self.geo_1)
        self.phys_1['pore.baz'] = 11
        self.phys_1['throat.baz'] = 11

        self.phys_2 = op.physics.GenericPhysics(network=self.net,
                                                phase=self.phase_1,
                                                geometry=self.geo_2)
        self.phys_2['pore.baz'] = 12
        self.phys_2['throat.baz'] = 12

        self.phys_3 = op.physics.GenericPhysics(network=self.net,
                                                phase=self.phase_2,
                                                geometry=self.geo_1)
        self.phys_3['pore.baz'] = 21
        self.phys_3['throat.baz'] = 21

        self.phys_4 = op.physics.GenericPhysics(network=self.net,
                                                phase=self.phase_2,
                                                geometry=self.geo_2)
        self.phys_4['pore.baz'] = 22
        self.phys_4['throat.baz'] = 22

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_save(self, tmpdir):
        fname = tmpdir.join(self.net.project.name)
        len_before = len(tmpdir.listdir())
        op.io.CSV.export_data(network=self.net, phases=self.phase_1, filename=fname)
        assert len(tmpdir.listdir()) == (len_before + 1)
        os.remove(fname.dirpath().join(self.net.project.name + '.csv'))

    def test_load_bad_filename(self, tmpdir):
        with pytest.raises(OSError):
            op.io.CSV.load(filename='')

    def test_load_categorized_by_object(self, tmpdir):
        fname = tmpdir.join(self.net.project.name)
        op.io.CSV.export_data(network=self.net, phases=self.phase_1, filename=fname)
        proj = op.io.CSV.import_data(filename=fname)
        os.remove(fname.dirpath().join(self.net.project.name + '.csv'))
        assert len(proj) == 2
        assert proj.network.name == self.net.name
        assert list(proj.phases().values())[0].name == self.phase_1.name


if __name__ == '__main__':
    # All the tests in this file can be run with 'playing' this file
    t = CSVTest()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=py.path.local())
