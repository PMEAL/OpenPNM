import numpy as np
from openpnm.io import GenericIO


class COMSOL(GenericIO):
    @classmethod
    def save(cls, network, phases=[], filename=''):
        project, network, phases = cls._parse_args(network=network,
                                                   phases=phases)
        network = network[0]
        f = open(filename+'.mphtxt', 'w')

        header(file=f, Nr=network.Nt, Nc=network.Np)

        cn = network['throat.conns']

        r_c = np.average(network['pore.coords'][cn], axis=1)
        r_w = network['throat.diameter']
        r_l = (network['throat.conduit_lengths.pore1'] +
               network['throat.conduit_lengths.throat'] +
               network['throat.conduit_lengths.pore2'])
        rectangles(file=f, centers=r_c, widths=r_w, lengths=r_l)

        c_c = network['pore.coords']
        c_r = network['pore.diameter']/2.0
        circles(file=f, centers=c_c, radii=c_r)
        print('OKK')

        f.close()


def rectangles(file, centers, widths, lengths):
    f = file
    for r in range(len(centers)):
        p1x = centers[r][0]-lengths[r]/2.0
        p1y = centers[r][1]-widths[r]/2.0
        p2x = centers[r][0]+lengths[r]/2.0
        p2y = centers[r][1]
        p3x = centers[r][0]+lengths[r]/2.0
        p3y = centers[r][1]+widths[r]/2.0
        p4x = centers[r][0]
        p4y = centers[r][1]+widths[r]/2.0

        f.write('# --------- rectangle object nbr ' +
                str(r+1)+' ----------'+' '+'\n')

        f.write('0 0 1'+' '+'\n')
        f.write('5 Geom2 # class'+' '+'\n')
        f.write('2 # version'+' '+'\n')
        f.write('2 # type'+' '+'\n')
        f.write('1 # voidsLabeled'+' '+'\n')
        f.write('1e-010 # gtol'+' '+'\n')
        f.write('0.0001 # resTol'+' '+'\n')
        f.write('4 # number of vertices'+' '+'\n')
        f.write('# Vertices'+' '+'\n')
        f.write('# X Y dom tol'+' '+'\n')

        f.write(str(p1x)+' '+str(p1y)+' '+'-1 NAN'+' '+'\n')
        f.write(str(p2x)+' '+str(p2y)+' '+'-1 NAN'+' '+'\n')
        f.write(str(p3x)+' '+str(p3y)+' '+'-1 NAN'+' '+'\n')
        f.write(str(p4x)+' '+str(p4y)+' '+'-1 NAN'+' '+'\n')

        f.write('4 # number of edges'+' '+'\n')
        f.write('# Edges'+' '+'\n')
        f.write('# vtx1 vtx2 s1 s2 up down mfd tol'+' '+'\n')
        f.write('2 1 0 1 0 1 1 NAN'+' '+'\n')
        f.write('3 2 0 1 0 1 2 NAN'+' '+'\n')
        f.write('4 3 0 1 0 1 3 NAN'+' '+'\n')
        f.write('1 4 0 1 0 1 4 NAN'+' '+'\n')
        f.write('4 # number of manifolds'+' '+'\n')
        f.write('# Manifolds'+' '+'\n')

        f.write('# Manifold #0'+' '+'\n')

        f.write('11 BezierCurve # class'+' '+'\n')
        f.write('0 0 # version'+' '+'\n')
        f.write('2 # sdim'+' '+'\n')
        f.write('0 2 1 # transformation'+' '+'\n')
        f.write('1 0 # degrees'+' '+'\n')
        f.write('2 # number of control points'+' '+'\n')
        f.write('# control point coords and weights'+' '+'\n')
        f.write(str(p2x)+' '+str(p2y)+' '+'1'+' '+'\n')
        f.write(str(p1x)+' '+str(p1y)+' '+'1'+' '+'\n')

        f.write('# Manifold #1'+' '+'\n')

        f.write('11 BezierCurve # class'+' '+'\n')
        f.write('0 0 # version'+' '+'\n')
        f.write('2 # sdim'+' '+'\n')
        f.write('0 2 1 # transformation'+' '+'\n')
        f.write('1 0 # degrees'+' '+'\n')
        f.write('2 # number of control points'+' '+'\n')
        f.write('# control point coords and weights'+' '+'\n')
        f.write(str(p3x)+' '+str(p3y)+' '+'1'+' '+'\n')
        f.write(str(p2x)+' '+str(p2y)+' '+'1'+' '+'\n')

        f.write('# Manifold #2'+' '+'\n')

        f.write('11 BezierCurve # class'+' '+'\n')
        f.write('0 0 # version'+' '+'\n')
        f.write('2 # sdim'+' '+'\n')
        f.write('0 2 1 # transformation'+' '+'\n')
        f.write('1 0 # degrees'+' '+'\n')
        f.write('2 # number of control points'+' '+'\n')
        f.write('# control point coords and weights'+' '+'\n')
        f.write(str(p4x)+' '+str(p4y)+' '+'1'+' '+'\n')
        f.write(str(p3x)+' '+str(p3y)+' '+'1'+' '+'\n')

        f.write('# Manifold #3'+' '+'\n')

        f.write('11 BezierCurve # class'+' '+'\n')
        f.write('0 0 # version'+' '+'\n')
        f.write('2 # sdim'+' '+'\n')
        f.write('0 2 1 # transformation'+' '+'\n')
        f.write('1 0 # degrees'+' '+'\n')
        f.write('2 # number of control points'+' '+'\n')
        f.write('# control point coords and weights'+' '+'\n')
        f.write(str(p1x)+' '+str(p1y)+' '+'1'+' '+'\n')
        f.write(str(p4x)+' '+str(p4y)+' '+'1'+' '+'\n')
        f.write('# Attributes'+' '+'\n')
        f.write('0 # nof attributes'+' '+'\n')

    return


def circles(file, centers, radii):
    f = file
    for c in range(len(centers)):
        p1x = centers[c][0]-radii[c]
        p1y = centers[c][1]
        p2x = centers[c][0]
        p2y = centers[c][1]-radii[c]
        p3x = centers[c][0]+radii[c]
        p3y = centers[c][1]
        p4x = centers[c][0]
        p4y = centers[c][1]+radii[c]

        f.write('# --------- Object circle ----------'+' '+'\n')
        f.write('# --------- circle object nbr ' +
                str(c+1)+' ----------'+' '+'\n')

        f.write('0 0 1'+' '+'\n')
        f.write('5 Geom2 # class'+' '+'\n')
        f.write('2 # version'+' '+'\n')
        f.write('2 # type'+' '+'\n')
        f.write('1 # voidsLabeled'+' '+'\n')
        f.write('1e-010 # gtol'+' '+'\n')
        f.write('0.0001 # resTol'+' '+'\n')
        f.write('4 # number of vertices'+' '+'\n')
        f.write('# Vertices'+' '+'\n')
        f.write('# X Y dom tol'+' '+'\n')

        # extreme left point of the circle
        f.write(str(p1x)+' '+str(p1y)+' '+'-1 NAN'+' '+'\n')
        # extreme bottom point of the circle
        f.write(str(p2x)+' '+str(p2y)+' '+'-1 NAN'+' '+'\n')
        # extreme right point of the circle
        f.write(str(p3x)+' '+str(p3y)+' '+'-1 NAN'+' '+'\n')
        # extreme top point of the circle
        f.write(str(p4x)+' '+str(p4y)+' '+'-1 NAN'+' '+'\n')

        f.write('4 # number of edges'+' '+'\n')
        f.write('# Edges'+' '+'\n')
        f.write('# vtx1 vtx2 s1 s2 up down mfd tol'+' '+'\n')
        f.write('2 1 0 1 0 1 1 NAN'+' '+'\n')
        f.write('3 2 0 1 0 1 2 NAN'+' '+'\n')
        f.write('4 3 0 1 0 1 3 NAN'+' '+'\n')
        f.write('1 4 0 1 0 1 4 NAN'+' '+'\n')
        f.write('4 # number of manifolds'+' '+'\n')
        f.write('# Manifolds'+' '+'\n')

        # bottom left quart of the circle
        f.write('# Manifold #0'+' '+'\n')

        f.write('11 BezierCurve # class'+' '+'\n')
        f.write('0 0 # version'+' '+'\n')
        f.write('2 # sdim'+' '+'\n')
        f.write('0 2 1 # transformation'+' '+'\n')
        f.write('2 0 # degrees'+' '+'\n')
        f.write('3 # number of control points'+' '+'\n')
        f.write('# control point coords and weights'+' '+'\n')
        # extreme bottom point of the circle
        f.write(str(p2x)+' '+str(p2y)+' '+'1'+' '+'\n')
        f.write(str(p1x)+' '+str(p2y)+' '+'0.70710678118654746'+' '+'\n')
        # extreme left point of the circle
        f.write(str(p1x)+' '+str(p1y)+' '+'1'+' '+'\n')

        # bottom right quart of the circle
        f.write('# Manifold #1'+' '+'\n')

        f.write('11 BezierCurve # class'+' '+'\n')
        f.write('0 0 # version'+' '+'\n')
        f.write('2 # sdim'+' '+'\n')
        f.write('0 2 1 # transformation'+' '+'\n')
        f.write('2 0 # degrees'+' '+'\n')
        f.write('3 # number of control points'+' '+'\n')
        f.write('# control point coords and weights'+' '+'\n')
        # extreme right point of the circle
        f.write(str(p3x)+' '+str(p3y)+' '+'1'+' '+'\n')
        f.write(str(p3x)+' '+str(p2y)+' '+'0.70710678118654746'+' '+'\n')
        # extreme bottom point of the circle
        f.write(str(p2x)+' '+str(p2y)+' '+'1'+' '+'\n')

        # top right quart of the circle
        f.write('# Manifold #2'+' '+'\n')

        f.write('11 BezierCurve # class'+' '+'\n')
        f.write('0 0 # version'+' '+'\n')
        f.write('2 # sdim'+' '+'\n')
        f.write('0 2 1 # transformation'+' '+'\n')
        f.write('2 0 # degrees'+' '+'\n')
        f.write('3 # number of control points'+' '+'\n')
        f.write('# control point coords and weights'+' '+'\n')
        # extreme top point of the circle
        f.write(str(p4x)+' '+str(p4y)+' '+'1'+' '+'\n')
        f.write(str(p3x)+' '+str(p4y)+' '+'0.70710678118654746'+' '+'\n')
        # extreme right point of the circle
        f.write(str(p3x)+' '+str(p3y)+' '+'1'+' '+'\n')

        # top left quart of the circle
        f.write('# Manifold #3'+' '+'\n')

        f.write('11 BezierCurve # class'+' '+'\n')
        f.write('0 0 # version'+' '+'\n')
        f.write('2 # sdim'+' '+'\n')
        f.write('0 2 1 # transformation'+' '+'\n')
        f.write('2 0 # degrees'+' '+'\n')
        f.write('3 # number of control points'+' '+'\n')
        f.write('# control point coords and weights'+' '+'\n')
        # extreme left point of the circle
        f.write(str(p1x)+' '+str(p1y)+' '+'1'+' '+'\n')
        f.write(str(p1x)+' '+str(p4y)+' '+'0.70710678118654746'+' '+'\n')
        # extreme top point of the circle
        f.write(str(p4x)+' '+str(p4y)+' '+'1'+' '+'\n')

        f.write('# Attributes'+' '+'\n')
        f.write('0 # nof attributes'+' '+'\n')

    return


def header(file, Nr, Nc):
    f = file

    f.write('# Geometry exported by OpenPNM'+'\n')

    f.write('# Major & minor version'+' '+'\n')
    f.write('0 1'+' '+'\n')
    f.write(str(Nc+Nr)+' '+'# number of tags'+' '+'\n')
    f.write('# Tags'+' '+'\n')

    for r in range(1, Nr+1):
        tag = 'sq'+str(int(r))
        tag = str(len(tag))+' '+tag
        f.write(tag+' '+'\n')

    for c in range(Nr+1, Nr+Nc+1):
        tag = 'c'+str(int(c))
        tag = str(len(tag))+' '+tag
        f.write(tag+' '+'\n')

    f.write(str(Nc+Nr)+' '+'# number of types'+' '+'\n')

    f.write('# Types'+' '+'\n')

    for i in range(Nc+Nr+1):
        f.write('3 obj'+' '+'\n')

    return
