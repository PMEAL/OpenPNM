import os as os
import numpy as np
from pathlib import Path
from openpnm.utils import logging, Project
from openpnm.network import GenericNetwork
from openpnm.io import GenericIO
from openpnm.topotools import trim
logger = logging.getLogger(__name__)


class MARock(GenericIO):
    r"""
    3DMA-Rock is a network extraction algorithm developed by Brent Lindquist
    and his group

    It uses Medial Axis thinning to find the skeleton of the pore space, then
    extracts geometrical features such as pore volume and throat
    cross-sectional area.

    [1] Lindquist, W. Brent, S. M. Lee, W. Oh, A. B. Venkatarangan, H. Shin,
    and M. Prodanovic. "3DMA-Rock: A software package for automated analysis
    of rock pore structure in 3-D computed microtomography images." SUNY Stony
    Brook (2005).
    """

    @classmethod
    def load(cls, *args, **kwargs):
        r"""
        This method is being deprecated.  Use ``import_data`` instead.
        """
        return cls.import_data(*args, **kwargs)

    @classmethod
    def import_data(cls, path, voxel_size=1, project=None):
        r"""
        Load data from a 3DMA-Rock extracted network.  This format consists of
        two files: 'rockname.np2th' and 'rockname.th2pn'.  They should be
        stored together in a folder which is referred to by the path argument.
        These files are binary and therefore not human readable.

        Parameters
        ----------
        path : string
            The location of the 'np2th' and 'th2np' files. This can be an
            absolute path or relative to the current working directory.

        network : OpenPNM Network Object
            If an Network object is recieved, this method will add new data to
            it but NOT overwrite anything that already exists.  This can be
            used to append data from different sources.

        voxel_size : scalar
            The resolution of the image on which 3DMA-Rock was run, in terms of
            the linear length of eac voxel. The default is 1.  This is used to
            scale the voxel counts to actual dimension. It is recommended that
            this value be in SI units [m] to work well with OpenPNM.

        project : OpenPNM Project object
            A GenericNetwork is created and added to the specified Project.
            If no Project is supplied then one will be created and returned.

        """

        net = {}
        path = Path(path)
        path = path.resolve()
        for file in os.listdir(path):
            if file.endswith(".np2th"):
                np2th_file = os.path.join(path, file)
            elif file.endswith(".th2np"):
                th2np_file = os.path.join(path, file)

        with open(np2th_file, mode='rb') as f:
            [Np, Nt] = np.fromfile(file=f, count=2, dtype='u4')
            net['pore.boundary_type'] = np.ndarray([Np, ], int)
            net['throat.conns'] = np.ones([Nt, 2], int)*(-1)
            net['pore.coordination'] = np.ndarray([Np, ], int)
            net['pore.ID_number'] = np.ndarray([Np, ], int)
            for i in range(0, Np):
                ID = np.fromfile(file=f, count=1, dtype='u4')
                net['pore.ID_number'][i] = ID
                net['pore.boundary_type'][i] = np.fromfile(file=f, count=1,
                                                           dtype='u1')
                z = np.fromfile(file=f, count=1, dtype='u4')[0]
                net['pore.coordination'][i] = z
                att_pores = np.fromfile(file=f, count=z, dtype='u4')
                att_throats = np.fromfile(file=f, count=z, dtype='u4')
                for j in range(0, len(att_throats)):
                    t = att_throats[j] - 1
                    p = att_pores[j] - 1
                    net['throat.conns'][t] = [i, p]
            net['throat.conns'] = np.sort(net['throat.conns'], axis=1)
            net['pore.volume'] = np.fromfile(file=f, count=Np, dtype='u4')
            nx = np.fromfile(file=f, count=1, dtype='u4')
            nxy = np.fromfile(file=f, count=1, dtype='u4')
            pos = np.fromfile(file=f, count=Np, dtype='u4')
            ny = nxy/nx
            ni = np.mod(pos, nx)
            nj = np.mod(np.floor(pos/nx), ny)
            nk = np.floor(np.floor(pos/nx)/ny)
            net['pore.coords'] = np.array([ni, nj, nk]).T

        with open(th2np_file, mode='rb') as f:
            Nt = np.fromfile(file=f, count=1, dtype='u4')[0]
            net['throat.area'] = np.ones([Nt, ], dtype=int)*(-1)
            for i in range(0, Nt):
                ID = np.fromfile(file=f, count=1, dtype='u4')
                net['throat.area'][i] = np.fromfile(file=f, count=1,
                                                    dtype='f4')
                # numvox = np.fromfile(file=f, count=1, dtype='u4')
                att_pores = np.fromfile(file=f, count=2, dtype='u4')
            nx = np.fromfile(file=f, count=1, dtype='u4')
            nxy = np.fromfile(file=f, count=1, dtype='u4')
            pos = np.fromfile(file=f, count=Nt, dtype='u4')
            ny = nxy/nx
            ni = np.mod(pos, nx)
            nj = np.mod(np.floor(pos/nx), ny)
            nk = np.floor(np.floor(pos/nx)/ny)
            net['throat.coords'] = np.array([ni, nj, nk]).T
            net['pore.internal'] = net['pore.boundary_type'] == 0

        # Convert voxel area and volume to actual dimensions
        net['throat.area'] = (voxel_size**2)*net['throat.area']
        net['pore.volume'] = (voxel_size**3)*net['pore.volume']

        if project is None:
            project = Project(name=path)
        network = GenericNetwork(project=project)
        network = cls._update_network(network=network, net=net)

        # Trim headless throats before returning
        ind = np.where(network['throat.conns'][:, 0] == -1)[0]
        trim(network=network, throats=ind)

        return project
