import pytest
import numpy as np
import openpnm as op
import matplotlib.pyplot as plt
from numpy.testing import assert_allclose


class PlotToolsTest:

    def setup_class(self):
        self.ws = op.Workspace()

    def test_plot_tutorial(self):
        pn = op.network.Cubic(shape=[4, 4, 1])
        # This runs locally, but fails on the CI due to a missing argument 's'
        # coming from within the networkx function, not ours.
        # op.visualization.plot_tutorial(pn)
        # plt.close()

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
            m = op.visualization.plot_networkx(pn, ax=ax)
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
            m = op.visualization.plot_networkx(pn, ax=ax)
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
            m = op.visualization.plot_networkx(pn, ax=ax)
            x_plot, y_plot = np.array(m.get_offsets()).T
            np.testing.assert_allclose(x_plot, x)
            plt.close()

    def test_plot_networkx_3d(self):
        pn = op.network.Cubic(shape=[5, 8, 3])
        with pytest.raises(Exception):
            op.visualization.plot_networkx(pn)

    def test_plot_connections_color_by(self):
        pn = op.network.Cubic(shape=[5, 5, 1])
        np.random.seed(10)
        Ts = np.array([0, 4, 6, 18])
        _ = op.visualization.plot_connections(
            pn,
            throats=Ts,
            color_by=np.random.rand(pn.Nt),
        )

    def test_plot_coordinates_color_by(self):
        pn = op.network.Cubic(shape=[5, 5, 1])
        np.random.seed(10)
        _ = op.visualization.plot_coordinates(
            pn,
            color_by=np.random.rand(pn.Np),
        )

    def test_plot_connections_and_coordinates_color_by_3D(self):
        pn = op.network.Cubic(shape=[5, 5, 5])
        np.random.seed(10)
        ax = op.visualization.plot_connections(
            pn,
            color_by=np.random.rand(pn.Nt),
        )
        ax = op.visualization.plot_coordinates(
            pn,
            color_by=np.random.rand(pn.Np),
            ax=ax,
        )

    def test_plot_connections_and_coordinates_color(self):
        pn = op.network.Cubic(shape=[5, 5, 1])
        ax = op.visualization.plot_connections(
            network=pn,
            color='b',
        )
        ax = op.visualization.plot_coordinates(
            network=pn,
            color='b',
            ax=ax,
        )

    def test_3d(self):
        pn = op.network.Cubic(shape=[10, 10, 3])
        pn['pore.internal'] = True
        pn.add_boundary_pores()
        Ps = pn.pores('internal')  # find internal pores
        fig, ax = plt.subplots()  # create empty figure
        ax = op.visualization.plot_coordinates(
            network=pn,
            pores=Ps,
            color='g',
            ax=ax,
        )
        Ps = pn.pores('*boundary')  # find boundary pores
        ax = op.visualization.plot_coordinates(
            network=pn,
            pores=Ps,
            color='r',
            ax=ax,
        )


if __name__ == '__main__':

    t = PlotToolsTest()
    self = t
    t.setup_class()
    for item in t.__dir__():
        if item.startswith('test'):
            print(f'Running test: {item}')
            t.__getattribute__(item)()
