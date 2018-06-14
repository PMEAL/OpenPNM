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
        shape = np.array(shape)
        if np.any(shape < 2):
            raise Exception('Bravais lattice networks must have at least 2 '
                            'pores in all directions')
        if mode == 'bcc':
            shape += 1
            net = CubicDual(shape=shape, spacing=spacing)
            self.update(net)
            topotools.trim(network=self, pores=net.pores(['left', 'right',
                                                          'front', 'back',
                                                          'top', 'bottom']))
            self['pore.coords'] -= 0.5
            # Deal with labels
            Ps1 = self['pore.secondary']
            self.clear(mode='labels')
            self['pore.corner_sites'] = Ps1
            self['pore.body_sites'] = ~Ps1
            Ts = self.find_neighbor_throats(pores=self.pores('body_sites'),
                                            mode='exclusive_or')
            self['throat.corner_to_body'] = False
            self['throat.corner_to_body'][Ts] = True
            Ts = self.find_neighbor_throats(pores=self.pores('corner_sites'),
                                            mode='intersection')
            self['throat.corner_to_corner'] = False
            self['throat.corner_to_corner'][Ts] = True
            Ts = self.find_neighbor_throats(pores=self.pores('body_sites'),
                                            mode='intersection')
            self['throat.body_to_body'] = False
            self['throat.body_to_body'][Ts] = True

        elif mode == 'fcc':
            shape = np.array(shape)
            # Create base cubic network of corner sites
            net1 = Cubic(shape=shape)
            net1['pore.corner_sites'] = True
            # Create 3 networks to become face sites
            net2 = Cubic(shape=shape - [1, 1, 0])
            net3 = Cubic(shape=shape - [1, 0, 1])
            net4 = Cubic(shape=shape - [0, 1, 1])
            net2['pore.coords'] += np.array([0.5, 0.5, 0])
            net3['pore.coords'] += np.array([0.5, 0, 0.5])
            net4['pore.coords'] += np.array([0, 0.5, 0.5])
            net2['pore.face_sites'] = True
            net3['pore.face_sites'] = True
            net4['pore.face_sites'] = True
            # Remove throats from net2 (trim doesn't work when removing ALL)
            for n in [net2, net3, net4]:
                n.clear(element='throat', mode='all')
                n.update({'throat.all': np.array([], dtype=bool)})
                n.update({'throat.conns': np.ndarray([0, 2], dtype=bool)})
            # Join networks 2, 3 and 4 into one with all face sites
            topotools.stitch(net2, net3, net2.Ps, net3.Ps,
                             len_min=0.70, len_max=0.75)
            topotools.stitch(net2, net4, net2.Ps, net4.Ps,
                             len_min=0.70, len_max=0.75)
            # Join face sites network with the corner sites network
            topotools.stitch(net1, net2, net1.Ps, net2.Ps,
                             len_min=0.70, len_max=0.75)
            self.update(net1)
            # Deal with labels
            Ps1 = self['pore.corner_sites']
            Ps2 = self['pore.face_sites']
            self.clear(mode='labels')
            self['pore.corner_sites'] = Ps1
            self['pore.face_sites'] = Ps2
            Ts = self.find_neighbor_throats(pores=self.pores('corner_sites'),
                                            mode='intersection')
            self['throat.corner_to_corner'] = False
            self['throat.corner_to_corner'][Ts] = True
            Ts = self.find_neighbor_throats(pores=self.pores('face_sites'))
            self['throat.corner_to_face'] = False
            self['throat.corner_to_face'][Ts] = True
        elif mode == 'hcp':
            raise Exception('hcp is not implemented yet')
        elif mode == 'sc':
            net = Cubic(shape=shape, spacing=spacing, **kwargs)
            self.update(net)
        else:
            raise Exception('Unrecognized lattice type: ' + mode)

        topotools.label_faces(self)
        self['pore.coords'] *= np.array(spacing)
