import numpy as np
from openpnm.io import GenericIO
from openpnm.utils import logging
logger = logging.getLogger(__name__)


class Salome(GenericIO):
    r"""
    Write a Salome (salome-platform.org) py script to generate a geometry.
    The exported py file should be loaded from Salome with "load script".
    """
    _header = '''
import numpy as np
def pnm_2_salome(cylinder_head, cylinder_tail, cylinder_r,
                 sphere_c, sphere_r, explicit=False):
    import numpy as np
    import math
    import GEOM
    import SALOMEDS
    import salome
    from salome.geom import geomBuilder
    salome.salome_init()
    geompy = geomBuilder.New()
    O = geompy.MakeVertex(0, 0, 0)
    OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
    OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
    OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
    s = len(sphere_r)
    c = len(cylinder_r)
    vertex_id = 0
    fuse_list = np.array([])
    for i in range(s):
        vertex_s = geompy.MakeVertex(sphere_c[i, 0], sphere_c[i, 1],
                                     sphere_c[i, 2])
        sphere = geompy.MakeSpherePntR(vertex_s, sphere_r[i])
        if explicit:
            geompy.addToStudy(vertex_s, "Vertex_{}".format(vertex_id))
            geompy.addToStudy(sphere, "Sphere_{}".format(i))
        vertex_id += 1
        fuse_list = np.append(fuse_list, sphere)
    for i in range(c):
        vertex_c1 = geompy.MakeVertex(cylinder_head[i, 0], cylinder_head[i, 1],
                                      cylinder_head[i, 2])
        vertex_c2 = geompy.MakeVertex(cylinder_tail[i, 0], cylinder_tail[i, 1],
                                      cylinder_tail[i, 2])
        line = geompy.MakeLineTwoPnt(vertex_c1, vertex_c2)
        length = abs(np.linalg.norm(cylinder_head[i] - cylinder_tail[i]))
        cylinder = geompy.MakeCylinder(vertex_c1, line, cylinder_r[i], length)
        if explicit:
            geompy.addToStudy(vertex_c1, "Vertex_{}".format(vertex_id))
        vertex_id += 1
        if explicit:
            geompy.addToStudy(vertex_c2, "Vertex_{}".format(vertex_id))
        vertex_id += 1
        if explicit:
            geompy.addToStudy(line, "Line_{}".format(i))
            geompy.addToStudy(cylinder, "Cylinder_{}".format(i))
        fuse_list = np.append(fuse_list, cylinder)
    fuse = geompy.MakeFuseList(list(fuse_list), True, True)
    geompy.addToStudy(fuse, 'network')
    geompy.addToStudy( O, 'O' )
    geompy.addToStudy( OX, 'OX' )
    geompy.addToStudy( OY, 'OY' )
    geompy.addToStudy( OZ, 'OZ' )
    if salome.sg.hasDesktop():
        salome.sg.updateObjBrowser()

    '''
    _footer = '''
pnm_2_salome(cylinder_head, cylinder_tail, cylinder_r,
             sphere_c, sphere_r, explicit=explicit)

    '''

    @classmethod
    def save(cls, *args, **kwargs):
        r"""
        This method is being deprecated.  Use ``export_data`` instead.
        """
        cls.export_data(*args, **kwargs)

    @classmethod
    def export_data(cls, network, phases=[], filename='', explicit=False):
        r"""
        Saves the network data and writes a Salome .py instruction file.

        Parameters
        ----------
        network : OpenPNM Network Object
            The network containing the desired data

        phases : list of OpenPNM Phase Objects (optional, default is none)
            A list of phase objects whose data are to be included

        Notes
        -----
        This method only saves the data, not any of the pore-scale models or
        other attributes.  To save an actual OpenPNM Project use the
        ``Workspace`` object.\n
        """
        project, network, phases = cls._parse_args(network=network,
                                                   phases=phases)
        net = network[0]
        explicit = explicit
        f = open(filename+'.py', 'w')
        f.write(cls._header+'\n')

        x = net['throat.endpoints.head']
        s = str(x.shape)
        f.write('cylinder_head = np.array([')
        np.savetxt(f, X=[x.flatten()], fmt='%.18e', delimiter=',', newline='')
        f.write('])\n')
        f.write('cylinder_head = np.reshape(cylinder_head, '+s+')\n')

        x = net['throat.endpoints.tail']
        s = str(x.shape)
        f.write('cylinder_tail = np.array([')
        np.savetxt(f, X=[x.flatten()], fmt='%.18e', delimiter=',', newline='')
        f.write('])\n')
        f.write('cylinder_tail = np.reshape(cylinder_tail, '+s+')\n')

        x = net['throat.diameter']/2
        s = str(x.shape)
        f.write('cylinder_r = np.array([')
        np.savetxt(f, X=[x.flatten()], fmt='%.18e', delimiter=',', newline='')
        f.write('])\n')
        f.write('cylinder_r = np.reshape(cylinder_r, '+s+')\n')

        x = net['pore.coords']
        s = str(x.shape)
        f.write('sphere_c = np.array([')
        np.savetxt(f, X=[x.flatten()], fmt='%.18e', delimiter=',', newline='')
        f.write('])\n')
        f.write('sphere_c = np.reshape(sphere_c, '+s+')\n')

        x = net['pore.diameter']/2
        s = str(x.shape)
        f.write('sphere_r = np.array([')
        np.savetxt(f, X=[x.flatten()], fmt='%.18e', delimiter=',', newline='')
        f.write('])\n')
        f.write('sphere_r = np.reshape(sphere_r, '+s+')\n')

        f.write('explicit = '+str(explicit))

        f.write(cls._footer)

        f.close()
