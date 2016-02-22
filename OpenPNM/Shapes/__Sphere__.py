# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 10:31:54 2016

@author: Tom
"""
import numpy as np
from transforms3d import _gohlketransforms as tr
from scipy.spatial import Delaunay

from OpenPNM.Base import logging
logger = logging.getLogger(__name__)


class Sphere():
    r"""
    Simple Cube Class
    """
    def __init__(self, rad=1.0, coord=[0,0,0],
                 num_theta=12.0, num_phi=12.0, **kwargs):
            r"""
            """
            self._rad = rad
            self._coord = np.asarray(coord)
            self._num_theta = num_theta
            self._num_phi = num_phi
    
    def stretch(self, direction='all', value=[0, 0, 0]):
        r"""
        Stretch the cube
        """
        value = np.asarray(value)
        if value <= 0.0 or value is None:
            logger.error("Please enter a positive non-zero value")
        if direction == 'x':
            self._dx *= value
        elif direction == 'y':
            self._dy *= value
        elif direction == 'z':
            self._dz *= value
        elif direction == 'all':
            if np.shape(value) == 3:
                self._dx *= value[0]
                self._dy *= value[1]
                self._dz *= value[2]
            else:
                self._dx *= value
                self._dy *= value
                self._dz *= value
        else:
            logger.error("Valid directions for stretch are 'x', 'y', 'z' or 'all'")
    
    def translate(self, direction='x', value=[0, 0, 0]):
        r"""
        Translate the cube
        """
        if direction == 'x':
            self._coord[0]+= value
        elif direction == 'y':
            self._coord[1] += value
        elif direction == 'z':
            self._coord[2] += value
        else:
            logger.error("Valid directions for translation are 'x', 'y', 'z'")
            
    def rotate(self, direction='all', value=[0, 0, 0]):
        r"""
        Rotate the cube
        """
        if direction == 'x':
            self._angle_x += value
        elif direction == 'y':
            self._angle_y += value
        elif direction == 'z':
            self._angle_z += value
        elif direction == 'all':
            if len(value) == 3:
                self._angle_x += value[0]
                self._angle_y += value[1]
                self._angle_z += value[2]
            else:
                logger.error("When direction is 'all' the value must be 3D")
    
    def compute_vertices(self):
        r"""
        Compute the vertices from the coordinates, dx, dy, dz and angles
        """
        # Firstly generate the 8 corner vertices of the cube as if they are centred about the origin
        verts = []
        for theta in np.linspace(-np.pi, np.pi, self._num_theta+1):
            for phi in np.linspace(-np.pi, np.pi, self._num_phi+1):
                verts.append([self._rad*np.cos(theta)*np.cos(phi),
                              self._rad*np.sin(theta)*np.cos(phi),
                              self._rad*np.sin(phi)])
        verts = np.asarray(verts)
                 
        # Now apply rotation about each principal axis
        #rotation_angles = [self._angle_x, self._angle_y, self._angle_z]
        #axes = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        #for a, x in zip(rotation_angles, axes):
        #    M = tr.rotation_matrix(a, x)
        #    verts = np.dot(verts, M[:3, :3].T)
        # Finally Translate
        verts += self._coord
        # Update vertices
        self._verts = verts
        #self._faces = verts[faces]

    def plot(self, fig=None):
        r"""
        Plot the cube
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection
        if fig is None:
            fig = plt.figure()
        ax = fig.gca(projection='3d')
        if 1 == 1:
            tri = Delaunay(self._verts, qhull_options='QJ Pp')
            ordered_verts = self._verts[tri.convex_hull]
            items = Poly3DCollection(ordered_verts, linewidths=1,
                                           alpha=0.2, zsort='min')
            trans = 1.0
            face_colours = [(1, 0, 0, trans),
                            (1, 0, 0, trans),
                            (0, 1, 0, trans),
                            (0, 1, 0, trans),
                            (0, 0, 1, trans),
                            (0, 0, 1, trans)]
            items.set_facecolor(face_colours)
            ax.add_collection(items)
        reset = False
        if reset:
            xmin,  ymin, zmin = np.min(self._verts, axis=0) - np.array([1,1,1])
            xmax,  ymax, zmax = np.max(self._verts, axis=0) + np.array([1,1,1])
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.set_zlim(zmin, zmax)
        else:
            lim = 5
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.set_zlim(-lim, lim)
        ax.scatter(self._coord[0], self._coord[1], self._coord[2], c='y')
        ax.scatter(self._verts[:, 0], self._verts[:, 1], self._verts[:, 2],c='r')

        ax.ticklabel_format(style='sci', scilimits=(0, 0))
        return fig