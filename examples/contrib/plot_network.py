import numpy as np
import sys
import vispy
import openpnm as op
from vispy import app, scene
from PyQt5 import QtWidgets
from matplotlib import cm


def create_pore_colors_from_array(a, cmap = 'viridis'):
    colormap = cm.get_cmap(cmap)
    return colormap(a/a.max())

def create_throat_colors_from_array(a, cmap = 'viridis'):
    colormap = cm.get_cmap(cmap)
    return np.repeat(colormap(a/a.max()),2,axis=0)


def create_visual_pnm(pn):
    canvas = scene.SceneCanvas(keys='interactive', show=True)
    view = canvas.central_widget.add_view()
    view.camera = 'turntable'
    view.camera.fov = 30
    view.camera.distance = 3*np.max(pn['pore.coords'])

    pore_color = create_pore_colors_from_array(pn['pore.diameter'], cmap = 'viridis')
    throat_color = create_throat_colors_from_array(pn['throat.diameter'], cmap = 'viridis')
    #plot spheres
    vis = scene.visuals.Markers(
        pos=pn['pore.coords'],
        size=pn['pore.diameter'],
        antialias=1,
        face_color = pore_color,
        edge_width=0,
        scaling=True,
        spherical=True,
    )
    vis.parent = view.scene
    #plot axis
    vispy.scene.visuals.XYZAxis(parent=view.scene)
    #set camera center
    view.camera.center = np.array((pn['pore.coords'][:,0].max()/2, pn['pore.coords'][:,1].max()/2, pn['pore.coords'][:,2].max()/2))
    #data preparation 
    lines = np.zeros((len(pn['throat.conns']), 2,3))
    for i in range(len(pn['throat.conns'])):
        pair = pn['throat.conns'][i]
        line = np.array([[pn['pore.coords'][pair[0]], pn['pore.coords'][pair[1]]]])
        lines[i,:,:] = line
    #plot throats
    vis2 = scene.visuals.Line(lines, width = 2, color = throat_color, connect='segments', antialias=True,)
    vis2.parent = view.scene
    

if __name__ == '__main__':
    import time
    app = QtWidgets.QApplication(sys.argv)
    pn = op.network.Cubic(shape=[30, 30, 30], spacing=[1, 1, 1], connectivity=6)
    geom = op.geometry.GenericGeometry(network=pn, pores=pn.Ps, throats=pn.Ts)
    geom['pore.diameter'] = np.random.rand(pn.Np)  
    print(f'Total number of pores: {pn.Np}')
    print(f'Total number of throats: {pn.Nt}')
    P12 = pn['throat.conns']          
    D12 = geom['pore.diameter'][P12]  
    Dt = np.amin(D12, axis=1)
    geom['throat.diameter'] = Dt
    h = pn.check_network_health()
    #op.topotools.trim(network=pn, pores=h['trim_pores'])

    create_visual_pnm(pn)

    app.exec_() 
