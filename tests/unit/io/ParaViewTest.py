import sys
import os
from os.path import realpath
from pathlib import Path
import openpnm as op
from openpnm.io.ParaView import export_data, open_paraview


class ParaViewTest():

    def setup_class(self):
        self.path = os.path.dirname(os.path.abspath(sys.argv[0]))

    def test_export_data(self):
        pn = op.network.Cubic(shape=[30, 40])
        geo = op.geometry.StickAndBall(network=pn, pores=pn.Ps, throats=pn.Ts)
        water = op.phases.Water(network=pn)
        _ = op.physics.Standard(network=pn, phase=water, geometry=geo)
        op.io.VTK.save(pn, water, 'test.vtp')
        export_data(pn, filename='test.vtp')
        os.remove('test.pvsm')
        os.remove('test.vtp')

    def test_open_paraview(self):
        path = Path(realpath(__file__), '../../../fixtures/VTK-VTP/test.pvsm')
        open_paraview(filename=path)


if __name__ == "__main__":
    t = ParaViewTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith("test"):
            print(f"Running test: {item}")
            t.__getattribute__(item)()
