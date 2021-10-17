import py
import os
import numpy as np
import openpnm as op


class HDF5Test:

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

        self.net['pore.object'] = np.ones(self.net.Np, dtype=object)

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_to_hdf5(self, tmpdir):
        fname = tmpdir.join(self.net.project.name)
        f = op.io.HDF5.export_data(network=[self.net],
                               phases=[self.phase_1, self.phase_2],
                               filename=fname)
        assert list(f.keys()) == ['net_01', 'phase_01', 'phase_02']
        filename = f.filename
        f.close()
        os.remove(filename)

    def test_print(self, tmpdir):
        fname = tmpdir.join(self.net.project.name)
        f = op.io.HDF5.export_data(network=[self.net], filename=fname,
                               interleave=False)
        op.io.HDF5.print_levels(f)
        op.io.HDF5.print_flattened(f)
        f.close()
        os.remove(fname.dirpath().join(self.net.project.name + '.hdf'))


if __name__ == '__main__':
    # All the tests in this file can be run with 'playing' this file
    t = HDF5Test()
    self = t  # For interacting with the tests at the command line
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print('running test: '+item)
            try:
                t.__getattribute__(item)()
            except TypeError:
                t.__getattribute__(item)(tmpdir=py.path.local())
