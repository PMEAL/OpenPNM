from openpnm.geometry import GenericGeometry
import scipy as sp


class Continuum(GenericGeometry):
    r"""
    - only works with cubic network

    """

    def __init__(self, **kwargs):

        super().__init__(**kwargs)
        # find network
        if 'network' in kwargs.keys():
            network = kwargs.pop('network')
        else:
            project = kwargs.pop('project')
            network = project.network
        # find spacing
        sp_x, sp_y, sp_z = network.spacing
        # attach pore properties
        self['pore.size_x'] = sp_x
        self['pore.size_y'] = sp_y
        self['pore.size_z'] = sp_z
        self['pore.area_x'] = sp_y*sp_z
        self['pore.area_y'] = sp_x*sp_z
        self['pore.area_z'] = sp_x*sp_y
        self['pore.volume'] = sp_x*sp_y*sp_z
        # attach thraot properties
        self['throat.size_x'] = sp_x
        self['throat.size_y'] = sp_y
        self['throat.size_z'] = sp_z
        # Find orientation of each throat and create a label
        conns = network['throat.conns']
        coords = network['pore.coords']
        temp = coords[conns]
        temp = sp.absolute(temp[:, 0] - temp[:, 1])
        network['throat.dir_x'] = temp[:, 0] > 0
        network['throat.dir_y'] = temp[:, 1] > 0
        network['throat.dir_z'] = temp[:, 2] > 0
        # initialize throat properties
        self['throat.length'] = 0
        self['throat.width'] = 0
        self['throat.height'] = 0
        self['throat.area'] = 0
        self['throat.volume'] = 0
        # for throats oriented in x-direction
        Ts = network.throats('throat.dir_x')
        self['throat.length'][Ts] = sp_x*1e-12 # approximately zero length
        self['throat.width'][Ts] = sp_y
        self['throat.height'][Ts] = sp_z
        self['throat.area'][Ts] = sp_y*sp_z
        self['throat.volume'][Ts] = sp_x*sp_y*sp_z*1e-12
        Ts = network.throats('throat.dir_y')
        self['throat.length'][Ts] = sp_y*1e-12
        self['throat.width'][Ts] = sp_x
        self['throat.height'][Ts] = sp_z
        self['throat.area'][Ts] = sp_x*sp_z
        self['throat.volume'][Ts] = sp_x*sp_y*sp_z*1e-12
        Ts = network.throats('throat.dir_z')
        self['throat.length'][Ts] = sp_z*1e-12
        self['throat.width'][Ts] = sp_x
        self['throat.height'][Ts] = sp_y
        self['throat.area'][Ts] = sp_x*sp_y
        self['throat.volume'][Ts] = sp_x*sp_y*sp_z*1e-12