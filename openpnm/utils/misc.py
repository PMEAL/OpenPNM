import scipy as _sp
import time as _time
from collections import OrderedDict


class PrintableList(list):
    def __str__(self):
        horizontal_rule = '―' * 78
        lines = [horizontal_rule]
        self.sort()
        for i, item in enumerate(self):
            if '._' not in item:
                lines.append('{0}\t: {1}'.format(i + 1, item))
        lines.append(horizontal_rule)
        return '\n'.join(lines)

    def __repr__(self):
        self.sort()
        return super().__repr__()


class PrintableDict(OrderedDict):
    def __init__(self, *args, **kwargs):
        self._header = 'value'
        if 'header' in kwargs:
            self._header = kwargs.pop('header')
        super().__init__(*args, **kwargs)

    def __repr__(self):
        text = dict(self).__str__()
        return text

    def __str__(self):
        header = '―' * 78
        lines = [header]
        lines.append('{0:<35s} {1}'.format('key', self._header))
        lines.append(header)
        for item in list(self.keys()):
            if type(self[item]) == _sp.ndarray:
                lines.append('{0:<35s} {1}'.format(item, _sp.shape(self[item])))
            else:
                lines.append('{0:<35s} {1}'.format(item, self[item]))
        lines.append(header)
        return '\n'.join(lines)


class HealthDict(PrintableDict):
    r"""
    This class adds a 'health' check to a standard dictionary.  This check
    looks into the dict values, and considers empty lists as healthy and all
    else as unhealthy.  If one or more entries is 'unhealthy' the health method
    returns False.
    """
    def __init__(self, header='status', **kwargs):
        super().__init__(header=header, **kwargs)

    def _get_health(self):
        health = True
        for item in list(self.keys()):
            if self[item] != []:
                health = False
        return health

    health = property(fget=_get_health)


def tic():
    r'''
    Homemade version of matlab tic and toc function, tic starts or resets
    the clock, toc reports the time since the last call of tic.
    '''
    global _startTime_for_tictoc
    _startTime_for_tictoc = _time.time()


def toc(quiet=False):
    r'''
    Homemade version of matlab tic and toc function, tic starts or resets
    the clock, toc reports the time since the last call of tic.

    Parameters
    ----------
    quiet : Boolean
        If False (default) then a message is output to the console.  If True
        the message is not displayed and the elapsed time is returned.
    '''
    if '_startTime_for_tictoc' in globals():
        t = _time.time() - _startTime_for_tictoc
        if quiet is False:
            print('Elapsed time in seconds: ', t)
        else:
            return t
    else:
        print("Toc: start time not set")


def unique_list(input_list):
    r"""
    For a given list (of points) remove any duplicates
    """
    output_list = []
    if len(input_list) > 0:
        dim = _sp.shape(input_list)[1]
        for i in input_list:
            match = False
            for j in output_list:
                if dim == 3:
                    if i[0] == j[0] and i[1] == j[1] and i[2] == j[2]:
                        match = True
                elif dim == 2:
                    if i[0] == j[0] and i[1] == j[1]:
                        match = True
                elif dim == 1:
                    if i[0] == j[0]:
                        match = True
            if match is False:
                output_list.append(i)
    return output_list


def amalgamate_data(objs=[], delimiter='_'):
    r"""
    Returns a dictionary containing ALL pore data from all netowrk and/or
    phase objects received as arguments

    Parameters
    ----------
    obj : list of OpenPNM objects
        The network and Phase objects whose data should be amalgamated into a
        single dict

    delimiter : string
        The delimiter to place between the prop name and the object name.  For
        instance \'pore.air_molar_density\' or \'pore.air|molar_density'\.  The
        use of underscores can be problematic for reloading the data since they
        are also used in multiple word properties.  The default is '_' for
        backwards compatibility, but the '|' option is preferred.

    Returns
    -------
    A standard Python dict containing all the data from the supplied OpenPNM
    objects
    """
    if type(objs) is not list:
        objs = list(objs)
    data_amalgamated = {}
    dlim = delimiter
    exclusion_list = ['pore.centroid', 'pore.vertices', 'throat.centroid',
                      'throat.offset_vertices', 'throat.vertices', 'throat.normal',
                      'throat.perimeter', 'pore.vert_index', 'throat.vert_index']
    for item in objs:
        mro = [module.__name__ for module in item.__class__.__mro__]
        # If Network object, combine Geometry and Network keys
        if 'GenericNetwork' in mro:
            keys = []
            for key in list(item.keys()):
                keys.append(key)
            for geom in item._geometries:
                for key in list(geom.keys()):
                    if key not in keys:
                        keys.append(key)
        else:
            if 'GenericPhase' in mro:
                keys = []
                for key in list(item.keys()):
                    keys.append(key)
                for physics in item._physics:
                    for key in list(physics.keys()):
                        if key not in keys:
                            keys.append(key)
        keys.sort()
        for key in keys:
            if key not in exclusion_list:
                try:
                    if _sp.amax(item[key]) < _sp.inf:
                        element = key.split('.')[0]
                        propname = key.split('.')[1]
                        dict_name = element + '.' + item.name + dlim + propname
                        if key in ['pore.coords', 'throat.conns',
                                   'pore.all', 'throat.all']:
                            dict_name = key
                        data_amalgamated.update({dict_name: item[key]})
                except TypeError:
                    pass
    return data_amalgamated


def conduit_lengths(network, throats=None, mode='pore'):
    r"""
    Return the respective lengths of the conduit components defined by the throat
    conns P1 T P2
    mode = 'pore' - uses pore coordinates
    mode = 'centroid' uses pore and throat centroids
    """
    if throats is None:
        throats = network.throats()
    Ps = network['throat.conns']
    pdia = network['pore.diameter']

    if mode == 'centroid':
        try:
            pcentroids = network['pore.centroid']
            tcentroids = network['throat.centroid']
            if _sp.sum(_sp.isnan(pcentroids)) + _sp.sum(_sp.isnan(tcentroids)) > 0:
                mode = 'pore'
            else:
                plen1 = _sp.sqrt(_sp.sum(_sp.square(pcentroids[Ps[:, 0]] -
                                         tcentroids), 1))-network['throat.length']/2
                plen2 = _sp.sqrt(_sp.sum(_sp.square(pcentroids[Ps[:, 1]] -
                                         tcentroids), 1))-network['throat.length']/2
        except KeyError:
            mode = 'pore'
    if mode == 'pore':
        # Find half-lengths of each pore
        pcoords = network['pore.coords']
        # Find the pore-to-pore distance, minus the throat length
        lengths = _sp.sqrt(_sp.sum(_sp.square(pcoords[Ps[:, 0]] -
                                   pcoords[Ps[:, 1]]), 1)) - network['throat.length']
        lengths[lengths < 0.0] = 2e-9
        # Calculate the fraction of that distance from the first pore
        try:
            fractions = pdia[Ps[:, 0]]/(pdia[Ps[:, 0]] + pdia[Ps[:, 1]])
            # Don't allow zero lengths
#            fractions[fractions == 0.0] = 0.5
#            fractions[fractions == 1.0] = 0.5
        except:
            fractions = 0.5
        plen1 = lengths*fractions
        plen2 = lengths*(1-fractions)

    return _sp.vstack((plen1, network['throat.length'], plen2)).T[throats]
