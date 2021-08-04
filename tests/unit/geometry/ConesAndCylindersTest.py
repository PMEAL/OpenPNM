import openpnm as op
from openpnm.geometry import ConesAndCylinders
import openpnm.models.geometry as gmods


class SpheresAndCylindersTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.geo = ConesAndCylinders(
            network=self.net, pores=self.net.Ps, throats=self.net.Ts
        )

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_models(self):
        models = self.geo.models
        assert models["pore.volume"]["model"] == gmods.pore_volume.sphere
        mod = gmods.throat_length.cones_and_cylinders
        assert models["throat.length"]["model"] == mod
        assert models["throat.volume"]["model"] == gmods.throat_volume.cylinder
        mod = gmods.throat_cross_sectional_area.cylinder
        assert models["throat.cross_sectional_area"]["model"] == mod
        mod = gmods.hydraulic_size_factors.cones_and_cylinders
        assert models["throat.hydraulic_size_factors"]["model"] == mod
        mod = gmods.diffusive_size_factors.cones_and_cylinders
        assert models["throat.diffusive_size_factors"]["model"] == mod


if __name__ == '__main__':

    t = SpheresAndCylindersTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
