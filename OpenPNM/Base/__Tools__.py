"""
###############################################################################
Tools:  Useful classes for use throughout the project
###############################################################################
"""
import scipy as _sp
from collections import OrderedDict as _odict


class PrintableList(list):
    def __str__(self):
        horizontal_rule = '-' * 60
        lines = [horizontal_rule]
        self.sort()
        for i, item in enumerate(self):
            lines.append('{0}\t: {1}'.format(i + 1, item))
        lines.append(horizontal_rule)
        return '\n'.join(lines)


class ObjectContainer(dict):
    r"""
    This dict is a kludgy (but awesome!) fix for the current way object handles
    are stored on other objects.  Currently ``net.geometries()`` returns a list
    containing the names (i.e. strings) of all **Geometry** objects (OR if a
    name is given it returns an actual handle to the object).  This behavior
    must be maintained for backwards compatibility, but it's a pain since you
    can't actually iterate over the list.  For instance, it would be nice if
    you could do ``[item.regenerate() for item in pn.geometries]``; however,
    the ``geometries`` attribute is a callable method and it returns a list of
    strings, which is useless.  I've been wanting to replace this with a *dict*
    for a while, but have finally figured out how to do it.  By making this
    dict "callable" it can return the list of strings as it does now AND also
    still be a regular dictionary. Eventually, we could remove the callable
    aspect (i.e. in V2.0).
    """
    def __call__(self, name=None):
        if self == {}:
            return []
        if name is None:
            objs = [item for item in self.keys()]
            if self[objs[0]]._isa('network'):
                objs = [self[objs[0]]]
        else:
            if type(name) is not list:
                name = [name]
            objs = [self[item] for item in name]
            if objs[0]._isa('network'):
                objs = objs[0]
        return objs


class PrintableDict(_odict):
    def __init__(self, *args, **kwargs):
        self._header = 'value'
        if 'header' in kwargs:
            self._header = kwargs.pop('header')
        super().__init__(*args, **kwargs)

    def __repr__(self):
        text = dict(self).__str__()
        return text

    def __str__(self):
        horizontal_rule = '-' * 60
        lines = [horizontal_rule]
        lines.append('{0:<25s} {1}'.format('key', self._header))
        lines.append(horizontal_rule)
        for item in list(self.keys()):
            lines.append('{0:<25s} {1}'.format(item, self[item]))
        lines.append(horizontal_rule)
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


class SetLocations():

    @staticmethod
    def add(obj, element, locations):
        net = obj._net
        element = obj._parse_element(element, single=True)
        locations = obj._parse_locations(locations)

        # Adapt the method depending on type of object received
        if obj._isa('Physics'):
            boss = obj.parent_phase
            objs = boss._physics
        elif obj._isa('Geometry'):
            boss = obj._net
            objs = net._geometries
        else:
            raise Exception('Setting locations does not apply to received object')

        # Ensure locations are not already assigned to another object
        temp = _sp.zeros(net._count(element=element), dtype=int)
        for item in objs:
            inds = boss._get_indices(element=element, labels=item.name)
            temp[inds] += 1
        temp[locations] += 1  # Increment proposed locations
        if _sp.any(temp[locations] > 1):
            raise Exception('Some of the given '+element+' are already ' +
                            'assigned to an existing object')

        # Create new 'all' label for new size
        new_len = obj._count(element=element) + _sp.size(locations)
        obj.update({element+'.all': _sp.ones((new_len, ), dtype=bool)})

        # Set locations in Network dictionary
        inds_orig = net._get_indices(element=element, labels=obj.name)
        if element+'.'+obj.name not in net.keys():
            net[element+'.'+obj.name] = False
        net[element+'.'+obj.name][locations] = True
        inds_new = net._get_indices(element=element, labels=obj.name)

        # Set locations in Phase dictionary (if obj is a Physics object)
        if obj._isa('Physics'):
            phase = obj.parent_phase
            if element+'.'+obj.name not in phase.keys():
                phase[element+'.'+obj.name] = False
            phase[element+'.'+obj.name][locations] = True

        # Increase size of labels (add False at new locations)
        labels = obj.labels()
        labels.remove(element+'.all')
        for item in labels:
            if item.split('.')[0] == element:
                net[element+'.'+'blank'] = False
                net[element+'.'+'blank'][inds_orig] = obj[item]
                obj[item] = net[element+'.'+'blank'][inds_new]
        net.pop(element+'.'+'blank', None)

    @staticmethod
    def drop(obj, element, locations):
        net = obj._net
        element = obj._parse_element(element, single=True)
        locations = obj._parse_locations(locations)

        obj_inds = net._map(element=element,
                            locations=locations,
                            target=obj)
        keep = ~obj._tomask(locations=obj_inds, element=element)
        for item in list(obj.keys()):
            if item.split('.')[0] == element:
                temp = obj[item][keep]
                obj.update({item: temp})

        # Drop locations from Network dictionary
        net[element+'.'+obj.name][locations] = False

        # Drop locations from Phase dictionary (if obj is a Physics)
        if obj._isa('Physics'):
            phase = obj.parent_phase
            phase[element+'.'+obj.name][locations] = False
