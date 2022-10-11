import os
import sys
from os.path import realpath
from pathlib import Path
import openpnm as op
import pytest


@pytest.mark.skip(reason="Installing paraview takes forever!")
class ParaViewTest:

    def setup_class(self):
        self.path = os.path.dirname(os.path.abspath(sys.argv[0]))

    def test_export_data(self):
        pn = op.network.Cubic(shape=[30, 40])
        _ = op.phase.Water(network=pn)
        op.io.project_to_vtk(project=pn.project, filename='test.vtp')
        op.io.project_to_paraview(project=pn.project, filename='test.vtp')
        os.remove('test.pvsm')
        os.remove('test.vtp')

    def test_open_paraview(self):
        path = Path(realpath(__file__), '../../../fixtures/VTK-VTP/test.pvsm')
        op.io.ParaView.open_paraview(filename=path)


if __name__ == "__main__":
    t = ParaViewTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith("test"):
            print(f"Running test: {item}")
            t.__getattribute__(item)()
