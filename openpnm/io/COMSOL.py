import numpy as np
from openpnm.io import GenericIO


class COMSOL(GenericIO):
    r"""
    Writes a file containing pores and throats of the network in a format that
    can be opened in COMSOL.

    Notes
    -----
    - The exported files contain COMSOL geometry objects, not meshes.
    - This class exports in 2D only.

    """
    @classmethod
    def save(cls, *args, **kwargs):
        r"""
        This method is being deprecated.  Use ``export_data`` instead.
        """
        cls.export_data(*args, **kwargs)

    @classmethod
    def export_data(cls, network, phases=[], filename=''):
        r"""
        Saves the network and geometry data from the given objects into the
        specified file. This exports in 2D only where throats and pores have
        rectangular and circular shapes, respectively.

        Parameters
        ----------
        network : OpenPNM Network Object
            The network containing the desired data

        phases : list of OpenPNM Phase Objects (optional, default is none)

        Notes
        -----
        This method only saves the network and geometry data, not any of the
        pore-scale models or other attributes.  To save an actual OpenPNM
        Project use the ``Workspace`` object.

        """
        project, network, phases = cls._parse_args(network=network,
                                                   phases=phases)
        network = network[0]
        f = open(filename+'.mphtxt', 'w')

        header(file=f, Nr=network.Nt, Nc=network.Np)

        cn = network['throat.conns']

        p1 = network['pore.coords'][cn[:, 0]]
        p2 = network['pore.coords'][cn[:, 1]]

        # Compute the rotation angle of throats
        dif_x = p2[:, 0]-p1[:, 0]
        dif_y = p2[:, 1]-p1[:, 1]
        # Avoid division by 0
        m = np.array([dif_x_i != 0 for dif_x_i in dif_x])
        r = np.zeros((len(dif_x)))
        r[~m] = np.inf
        r[m] = dif_y[m]/dif_x[m]
        angles = np.arctan(r)

        r_w = network['throat.diameter']
        rectangles(file=f, pores1=p1, pores2=p2, alphas=angles, widths=r_w)

        c_c = network['pore.coords']
        c_r = network['pore.diameter']/2.0
        circles(file=f, centers=c_c, radii=c_r)

        f.close()


def header(file, Nr, Nc):
    f = file

    f.write('# Geometry exported by OpenPNM'+2*'\n')

    f.write('# Major & minor version'+'\n')
    f.write('0 1'+2*'\n')

    f.write(str(Nc+Nr)+' '+'# number of tags'+'\n')
    f.write('# Tags'+'\n')

    for r in range(1, Nr+1):
        tag = 'r'+str(int(r))
        tag = str(len(tag))+' '+tag
        f.write(tag+'\n')

    for c in range(1, Nc+1):
        tag = 'c'+str(int(c))
        tag = str(len(tag))+' '+tag
        f.write(tag+'\n')

    f.write('\n'+str(Nc+Nr)+' '+'# number of types'+'\n')
    f.write('# Types'+'\n')

    for i in range(Nc+Nr):
        f.write('3 obj'+'\n')

    f.write('\n')


def rectangles(file, pores1, pores2, alphas, widths):
    f = file

    p1x = pores1[:, 0] + (widths/2)*np.sin(alphas)
    p1y = pores1[:, 1] - (widths/2)*np.cos(alphas)
    p2x = pores2[:, 0] + (widths/2)*np.sin(alphas)
    p2y = pores2[:, 1] - (widths/2)*np.cos(alphas)
    p3x = pores2[:, 0] - (widths/2)*np.sin(alphas)
    p3y = pores2[:, 1] + (widths/2)*np.cos(alphas)
    p4x = pores1[:, 0] - (widths/2)*np.sin(alphas)
    p4y = pores1[:, 1] + (widths/2)*np.cos(alphas)

    for r in range(len(pores1)):
        f.write('# --------- rectangle nbr '+str(r+1)+' ---------'+2*'\n')

        f.write('0 0 1'+'\n')
        f.write('5 Geom2 # class'+'\n')
        f.write('2 # version'+'\n')
        f.write('2 # type'+'\n')
        f.write('1 # voidsLabeled'+'\n')
        f.write('1e-010 # gtol'+'\n')
        f.write('0.0001 # resTol'+2*'\n')

        f.write('4 # number of vertices'+'\n')
        f.write('# Vertices'+'\n')
        f.write('# X Y dom tol'+'\n')
        f.write(str(p1x[r])+' '+str(p1y[r])+' '+'-1 NAN'+'\n')
        f.write(str(p2x[r])+' '+str(p2y[r])+' '+'-1 NAN'+'\n')
        f.write(str(p3x[r])+' '+str(p3y[r])+' '+'-1 NAN'+'\n')
        f.write(str(p4x[r])+' '+str(p4y[r])+' '+'-1 NAN'+2*'\n')

        f.write('4 # number of edges'+'\n')
        f.write('# Edges'+'\n')
        f.write('# vtx1 vtx2 s1 s2 up down mfd tol'+'\n')
        f.write('2 1 0 1 0 1 1 NAN'+'\n')
        f.write('3 2 0 1 0 1 2 NAN'+'\n')
        f.write('4 3 0 1 0 1 3 NAN'+'\n')
        f.write('1 4 0 1 0 1 4 NAN'+2*'\n')

        f.write('4 # number of manifolds'+'\n')
        f.write('# Manifolds'+2*'\n')

        f.write('# Manifold #0'+'\n')
        f.write('11 BezierCurve # class'+'\n')
        f.write('0 0 # version'+'\n')
        f.write('2 # sdim'+'\n')
        f.write('0 2 1 # transformation'+'\n')
        f.write('1 0 # degrees'+'\n')
        f.write('2 # number of control points'+'\n')
        f.write('# control point coords and weights'+'\n')
        f.write(str(p2x[r])+' '+str(p2y[r])+' '+'1'+'\n')
        f.write(str(p1x[r])+' '+str(p1y[r])+' '+'1'+2*'\n')

        f.write('# Manifold #1'+'\n')
        f.write('11 BezierCurve # class'+'\n')
        f.write('0 0 # version'+'\n')
        f.write('2 # sdim'+'\n')
        f.write('0 2 1 # transformation'+'\n')
        f.write('1 0 # degrees'+'\n')
        f.write('2 # number of control points'+'\n')
        f.write('# control point coords and weights'+'\n')
        f.write(str(p3x[r])+' '+str(p3y[r])+' '+'1'+'\n')
        f.write(str(p2x[r])+' '+str(p2y[r])+' '+'1'+2*'\n')

        f.write('# Manifold #2'+'\n')
        f.write('11 BezierCurve # class'+'\n')
        f.write('0 0 # version'+'\n')
        f.write('2 # sdim'+'\n')
        f.write('0 2 1 # transformation'+'\n')
        f.write('1 0 # degrees'+'\n')
        f.write('2 # number of control points'+'\n')
        f.write('# control point coords and weights'+'\n')
        f.write(str(p4x[r])+' '+str(p4y[r])+' '+'1'+'\n')
        f.write(str(p3x[r])+' '+str(p3y[r])+' '+'1'+2*'\n')

        f.write('# Manifold #3'+'\n')
        f.write('11 BezierCurve # class'+'\n')
        f.write('0 0 # version'+'\n')
        f.write('2 # sdim'+'\n')
        f.write('0 2 1 # transformation'+'\n')
        f.write('1 0 # degrees'+'\n')
        f.write('2 # number of control points'+'\n')
        f.write('# control point coords and weights'+'\n')
        f.write(str(p1x[r])+' '+str(p1y[r])+' '+'1'+'\n')
        f.write(str(p4x[r])+' '+str(p4y[r])+' '+'1'+2*'\n')

        f.write('# Attributes'+'\n')
        f.write('0 # nof attributes'+2*'\n')


def circles(file, centers, radii):
    f = file

    p1x = centers[:, 0]-radii
    p1y = centers[:, 1]
    p2x = centers[:, 0]
    p2y = centers[:, 1]-radii
    p3x = centers[:, 0]+radii
    p3y = centers[:, 1]
    p4x = centers[:, 0]
    p4y = centers[:, 1]+radii

    for c in range(len(centers)):
        f.write('# --------- circle nbr '+str(c+1)+' ---------'+2*'\n')

        f.write('0 0 1'+'\n')
        f.write('5 Geom2 # class'+'\n')
        f.write('2 # version'+'\n')
        f.write('2 # type'+'\n')
        f.write('1 # voidsLabeled'+'\n')
        f.write('1e-010 # gtol'+'\n')
        f.write('0.0001 # resTol'+2*'\n')

        f.write('4 # number of vertices'+'\n')
        f.write('# Vertices'+'\n')
        f.write('# X Y dom tol'+'\n')
        # extreme left point of the circle
        f.write(str(p1x[c])+' '+str(p1y[c])+' '+'-1 NAN'+'\n')
        # extreme bottom point of the circle
        f.write(str(p2x[c])+' '+str(p2y[c])+' '+'-1 NAN'+'\n')
        # extreme right point of the circle
        f.write(str(p3x[c])+' '+str(p3y[c])+' '+'-1 NAN'+'\n')
        # extreme top point of the circle
        f.write(str(p4x[c])+' '+str(p4y[c])+' '+'-1 NAN'+2*'\n')

        f.write('4 # number of edges'+'\n')
        f.write('# Edges'+'\n')
        f.write('# vtx1 vtx2 s1 s2 up down mfd tol'+'\n')
        f.write('2 1 0 1 0 1 1 NAN'+'\n')
        f.write('3 2 0 1 0 1 2 NAN'+'\n')
        f.write('4 3 0 1 0 1 3 NAN'+'\n')
        f.write('1 4 0 1 0 1 4 NAN'+2*'\n')

        f.write('4 # number of manifolds'+'\n')
        f.write('# Manifolds'+2*'\n')

        # bottom left quart of the circle
        f.write('# Manifold #0'+'\n')
        f.write('11 BezierCurve # class'+'\n')
        f.write('0 0 # version'+'\n')
        f.write('2 # sdim'+'\n')
        f.write('0 2 1 # transformation'+'\n')
        f.write('2 0 # degrees'+'\n')
        f.write('3 # number of control points'+'\n')
        f.write('# control point coords and weights'+'\n')
        # extreme bottom point of the circle
        f.write(str(p2x[c])+' '+str(p2y[c])+' '+'1'+'\n')
        f.write(str(p1x[c])+' '+str(p2y[c])+' '+'0.70710678118654746'+'\n')
        # extreme left point of the circle
        f.write(str(p1x[c])+' '+str(p1y[c])+' '+'1'+2*'\n')

        # bottom right quart of the circle
        f.write('# Manifold #1'+'\n')
        f.write('11 BezierCurve # class'+'\n')
        f.write('0 0 # version'+'\n')
        f.write('2 # sdim'+'\n')
        f.write('0 2 1 # transformation'+'\n')
        f.write('2 0 # degrees'+'\n')
        f.write('3 # number of control points'+'\n')
        f.write('# control point coords and weights'+'\n')
        # extreme right point of the circle
        f.write(str(p3x[c])+' '+str(p3y[c])+' '+'1'+'\n')
        f.write(str(p3x[c])+' '+str(p2y[c])+' '+'0.70710678118654746'+'\n')
        # extreme bottom point of the circle
        f.write(str(p2x[c])+' '+str(p2y[c])+' '+'1'+2*'\n')

        # top right quart of the circle
        f.write('# Manifold #2'+'\n')
        f.write('11 BezierCurve # class'+'\n')
        f.write('0 0 # version'+'\n')
        f.write('2 # sdim'+'\n')
        f.write('0 2 1 # transformation'+'\n')
        f.write('2 0 # degrees'+'\n')
        f.write('3 # number of control points'+'\n')
        f.write('# control point coords and weights'+'\n')
        # extreme top point of the circle
        f.write(str(p4x[c])+' '+str(p4y[c])+' '+'1'+'\n')
        f.write(str(p3x[c])+' '+str(p4y[c])+' '+'0.70710678118654746'+'\n')
        # extreme right point of the circle
        f.write(str(p3x[c])+' '+str(p3y[c])+' '+'1'+2*'\n')

        # top left quart of the circle
        f.write('# Manifold #3'+'\n')
        f.write('11 BezierCurve # class'+'\n')
        f.write('0 0 # version'+'\n')
        f.write('2 # sdim'+'\n')
        f.write('0 2 1 # transformation'+'\n')
        f.write('2 0 # degrees'+'\n')
        f.write('3 # number of control points'+'\n')
        f.write('# control point coords and weights'+'\n')
        # extreme left point of the circle
        f.write(str(p1x[c])+' '+str(p1y[c])+' '+'1'+'\n')
        f.write(str(p1x[c])+' '+str(p4y[c])+' '+'0.70710678118654746'+'\n')
        # extreme top point of the circle
        f.write(str(p4x[c])+' '+str(p4y[c])+' '+'1'+2*'\n')

        f.write('# Attributes'+'\n')
        f.write('0 # nof attributes'+2*'\n')
