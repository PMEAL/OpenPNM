import openpnm as op
from openpnm.geometry import TrapezoidsAndRectangles
import openpnm.models.geometry as gmods


class SquaresAndRectanglesTest:

    def setup_class(self):
        self.net = op.network.Cubic(shape=[3, 3, 3])
        self.geo = TrapezoidsAndRectangles(
            network=self.net, pores=self.net.Ps, throats=self.net.Ts
        )

    def teardown_class(self):
        ws = op.Workspace()
        ws.clear()

    def test_models(self):
        models = self.geo.models
        assert models["pore.volume"]["model"] == gmods.pore_volume.circle
        mod = gmods.throat_length.trapezoids_and_rectangles
        assert models["throat.length"]["model"] == mod
        mod = gmods.throat_cross_sectional_area.rectangle
        assert models["throat.cross_sectional_area"]["model"] == mod
        assert models["throat.volume"]["model"] == gmods.throat_volume.rectangle
        mod = gmods.hydraulic_size_factors.trapezoids_and_rectangles
        assert models["throat.hydraulic_size_factors"]["model"] == mod
        mod = gmods.diffusive_size_factors.trapezoids_and_rectangles
        assert models["throat.diffusive_size_factors"]["model"] == mod


if __name__ == '__main__':

    t = SquaresAndRectanglesTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
