"""
===============================================================================
Bravais:
===============================================================================

"""
from openpnm.network import GenericNetwork, Cubic, CubicDual
from openpnm import topotools
from openpnm.core import logging
import numpy as np
logger = logging.getLogger(__name__)


class Bravais(GenericNetwork):
    r"""

    """
    def __init__(self, shape, mode, spacing=1, **kwargs):
        super().__init__(**kwargs)
        if mode == 'bcc':
            shape = np.array(shape) + 1
            if np.any(shape < 2):
                raise Exception('The bcc network must have at least 2 pores' +
                                'in all directions')
            net = CubicDual(shape=shape, spacing=spacing)
            self.update(net)
            topotools.trim(network=self, pores=net.pores(['left', 'right',
                                                          'front', 'back',
                                                          'top', 'bottom']))
            self.clear(mode='labels')
        elif mode == 'fcc':
            shape = np.array(shape)
            net1 = Cubic(shape=shape)
            net1['pore.corner_sites'] = True
            net2 = Cubic(shape=shape - [1, 1, 0])
            net2['pore.coords'] += np.array([0.5, 0.5, 0])
            net3 = Cubic(shape=shape - [1, 0, 1])
            net3['pore.coords'] += np.array([0.5, 0, 0.5])
            net4 = Cubic(shape=shape - [0, 1, 1])
            net4['pore.coords'] += np.array([0, 0.5, 0.5])
            topotools.extend(network=net2,
                             pore_coords=net3['pore.coords'])
            topotools.extend(network=net2,
                             pore_coords=net4['pore.coords'])
            net2['pore.face_sites'] = True
            topotools.stitch(net1, net2, net1.Ps, net2.Ps, len_min=.7,
                             len_max=.75)
            self.update(net1)
            Ts = self.find_neighbor_throats(pores=self.pores('face_sites'),
                                            mode='intersection')
            topotools.trim(network=self, throats=Ts)
        elif mode == 'hcp':
            pass
        elif mode == 'sc':
            net = CubicDual(shape=shape, spacing=spacing)
            self.update(net)
        else:
            raise Exception('Unrecognized lattice type: ' + mode)

        topotools.label_faces(self)
        self['pore.coords'] *= np.array(spacing)
