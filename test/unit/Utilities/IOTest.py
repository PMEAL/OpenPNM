import OpenPNM as op
import scipy as sp
import os
import OpenPNM.Utilities.IO as io


class IOTest:
    def setup_class(self):
        self.net = op.Network.Cubic(shape=[3, 3, 3])
        self.geom = op.Geometry.Cube_and_Cuboid(network=self.net,
                                                pores=self.net.Ps,
                                                throats=self.net.Ts)
        self.phase = op.Phases.Air(network=self.net)
        self.physics = op.Physics.Standard(network=self.net,
                                           phase=self.phase,
                                           pores=self.net.Ps,
                                           throats=self.net.Ts)

    def test_save_vtk_no_phases(self):
        fname = os.path.join(TEMP_DIR, 'test_save_vtk_1')
        io.VTK.save(network=self.net, filename=fname)

    def test_save_vtk_w_phases(self):
        fname = os.path.join(TEMP_DIR, 'test_save_vtk_2')
        io.VTK.save(network=self.net, filename=fname, phases=self.phase)

    def test_save_csv_no_phases(self):
        fname = os.path.join(TEMP_DIR, 'test_save_csv_1')
        io.CSV.save(network=self.net, filename=fname)

    def test_save_csv_w_phases(self):
        fname = os.path.join(TEMP_DIR, 'test_save_csv_2')
        io.CSV.save(network=self.net, filename=fname, phases=self.phase)

    def test_save_mat_no_phases(self):
        fname = os.path.join(TEMP_DIR, 'test_save_mat_1')
        io.MAT.save(network=self.net, filename=fname)

    def test_save_mat_w_phases(self):
        fname = os.path.join(TEMP_DIR, 'test_save_mat_2')
        io.MAT.save(network=self.net, filename=fname, phases=self.phase)
