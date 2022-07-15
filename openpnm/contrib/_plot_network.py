import numpy as np
import openpnm as op
from matplotlib import cm


__all__ = [
    'plot_3d',
]


def create_pore_colors_from_array(a, cmap='viridis'):
    colormap = cm.get_cmap(cmap)
    return colormap(a/a.max())


def create_throat_colors_from_array(a, cmap='viridis'):
    colormap = cm.get_cmap(cmap)
    return np.repeat(colormap(a/a.max()), 2, axis=0)


def plot_3d(
    network,
    pore_color=None,
    pore_size=None,
    throat_color=None,
    throat_size=None,
    bgcolor='grey',
):
    r"""

    Parameters
    ----------
    network
    """
    try:
        from vispy import scene
    except ModuleNotFoundError:
        raise Exception("vispy must be installed to use this function")
    canvas = scene.SceneCanvas(keys='interactive', show=True, bgcolor=bgcolor)
    view = canvas.central_widget.add_view()
    view.camera = 'turntable'
    view.camera.fov = 30
    view.camera.distance = 3*np.max(network['pore.coords'])

    if pore_color is None:
        pore_color = create_pore_colors_from_array(network['pore.diameter'],
                                                   cmap='viridis')
    else:
        pore_color = create_pore_colors_from_array(pore_color,
                                                   cmap='viridis')
    if throat_color is None:
        throat_color = create_throat_colors_from_array(network['throat.diameter'],
                                                       cmap='viridis')
    else:
        throat_color = create_throat_colors_from_array(throat_color,
                                                       cmap='viridis')
    if pore_size is None:
        pore_size = network['pore.diameter']
    if throat_size is None:
        throat_size = 2
    else:
        throat_size = np.max(throat_size)  # Arrays not supported here
    # plot spheres
    vis = scene.visuals.Markers(
        pos=network['pore.coords'],
        size=pore_size,
        antialias=0,
        face_color=pore_color,
        edge_width=0,
        scaling=True,
        spherical=True,
    )
    vis.parent = view.scene
    # plot axis
    # vispy.scene.visuals.XYZAxis(parent=view.scene)
    # set camera center
    view.camera.center = np.array((network['pore.coords'][:, 0].max()/2,
                                   network['pore.coords'][:, 1].max()/2,
                                   network['pore.coords'][:, 2].max()/2))
    # data preparation
    lines = np.zeros((len(network['throat.conns']), 2, 3))
    for i in range(len(network['throat.conns'])):
        pair = network['throat.conns'][i]
        line = np.array([[network['pore.coords'][pair[0]],
                          network['pore.coords'][pair[1]]]])
        lines[i, :, :] = line
    # plot throats
    vis2 = scene.visuals.Line(lines,
                              width=throat_size,
                              color=throat_color,
                              connect='segments',
                              antialias=True,)
    vis2.parent = view.scene


if __name__ == '__main__':
    pn = op.network.Cubic(shape=[30, 30, 30], spacing=1)
    pn['pore.diameter'] = np.random.rand(pn.Np)
    print(f'Total number of pores: {pn.Np}')
    print(f'Total number of throats: {pn.Nt}')
    P12 = pn['throat.conns']
    D12 = pn['pore.diameter'][P12]
    Dt = np.amin(D12, axis=1)
    pn['throat.diameter'] = Dt
    plot_3d(network=pn,
            pore_color=pn['pore.diameter'],
            throat_color=pn['throat.diameter'],
            throat_size=10*pn['throat.diameter'],
            bgcolor='grey')
