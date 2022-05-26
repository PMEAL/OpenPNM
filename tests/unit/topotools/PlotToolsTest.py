import pytest
import numpy as np
import openpnm as op
import matplotlib.pyplot as plt
from matplotlib import cm
from openpnm import topotools
from numpy.testing import assert_allclose


class PlotToolsTest:

    def setup_class(self):
        self.ws = op.Workspace()

    def test_plot_tutorial(self):
        pn = op.network.Cubic(shape=[4, 4, 1])
        topotools.plot_tutorial(pn)
        plt.close()

    def test_plot_networkx_var_spacing(self):
        for i in range(3):
            shape = np.ones(3, dtype=int)
            shape[np.arange(3) != i] = [5, 8]
            spacing = np.ones(3, dtype=float)
            spacing[np.arange(3) != i] = [0.01, 0.6]
            pn = op.network.Cubic(shape=shape)
            dims = op.topotools.dimensionality(pn)
            x, y = pn["pore.coords"].T[dims]
            fig, ax = plt.subplots()
            m = op.topotools.plot_networkx(pn, ax=ax)
            x_plot, y_plot = np.array(m.get_offsets()).T
            np.testing.assert_allclose(x_plot, x)
            np.testing.assert_allclose(y_plot, y)
            plt.close()

    def test_plot_networkx(self):
        # 2D networks in XY, YZ, XZ planes
        for i in range(3):
            shape = np.ones(3, dtype=int)
            shape[np.arange(3) != i] = [5, 8]
            pn = op.network.Cubic(shape=shape)
            x, y = pn["pore.coords"].T[op.topotools.dimensionality(pn)]
            fig, ax = plt.subplots()
            m = op.topotools.plot_networkx(pn, ax=ax)
            x_plot, y_plot = np.array(m.get_offsets()).T
            np.testing.assert_allclose(x_plot, x)
            np.testing.assert_allclose(y_plot, y)
            plt.close()
        # 1D networks in XY, YZ, XZ planes
        for i in range(3):
            shape = np.ones(3, dtype=int)
            shape[np.arange(3) == i] = [5]
            pn = op.network.Cubic(shape=shape)
            x, = pn["pore.coords"].T[op.topotools.dimensionality(pn)]
            fig, ax = plt.subplots()
            m = op.topotools.plot_networkx(pn, ax=ax)
            x_plot, y_plot = np.array(m.get_offsets()).T
            np.testing.assert_allclose(x_plot, x)
            plt.close()

    def test_plot_networkx_3d(self):
        pn = op.network.Cubic(shape=[5, 8, 3])
        with pytest.raises(Exception):
            op.topotools.plot_networkx(pn)

    def test_generate_voxel_image(self):
        pn = op.network.Cubic(shape=[5, 5, 1])
        pn.add_model_collection(
            op.models.collections.geometry.spheres_and_cylinders()
        )
        pn.regenerate_models()
        im = op.topotools.generate_voxel_image(network=pn,
                                               pore_shape='sphere',
                                               throat_shape='cylinder',
                                               max_dim=500)
        assert im.shape[0] == 500

    def test_plot_connections_color_by(self):
        pn = op.network.Cubic(shape=[5, 5, 1])
        np.random.seed(10)
        pn.add_model_collection(
            op.models.collections.geometry.spheres_and_cylinders()
        )
        pn.regenerate_models()
        Ts = np.array([0, 4, 6, 18])
        im = op.topotools.plot_connections(pn, throats=Ts,
                                           color_by=pn['throat.diameter'])
        colors_im = im.get_color()
        color_by = pn['throat.diameter'][Ts]
        color_calc = cm.get_cmap(name='jet')(color_by / color_by.max())
        color_calc[:, 3] = 1.0
        assert_allclose(color_calc, colors_im, rtol=1e-5)

    def test_plot_coordinates_color_by(self):
        pn = op.network.Cubic(shape=[5, 5, 1])
        np.random.seed(10)
        pn.add_model_collection(
            op.models.collections.geometry.spheres_and_cylinders()
        )
        pn.regenerate_models()
        Ps = np.array([0, 4, 6, 18])
        im = op.topotools.plot_coordinates(pn, pores=Ps,
                                           color_by=pn['pore.diameter'])
        colors_im = im.get_edgecolors()
        color_by = pn['pore.diameter'][Ps]
        color_calc = cm.get_cmap(name='jet')(color_by / color_by.max())
        color_calc[:, 3] = 1.0
        assert_allclose(color_calc, colors_im, rtol=1e-5)


if __name__ == '__main__':

    t = PlotToolsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
