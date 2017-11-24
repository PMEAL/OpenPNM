"""
###############################################################################
tools:  Useful classes for use throughout the project
###############################################################################
"""
from collections import OrderedDict
import scipy as sp

class PrintableList(list):
    def __str__(self):
        horizontal_rule = '―' * 60
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
        horizontal_rule = '―' * 60
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


def add_indices(obj, element, indices):
    sim = obj.simulation
    element = obj._parse_element(element, single=True)
    indices = obj._parse_indices(indices)

    # Adapt the method depending on type of object received
    if 'GenericPhysics' in obj.mro():
        boss = sim.find_phase(obj)
        objs = sim.physics.values()
    elif 'GenericGeometry' in obj.mro():
        boss = sim.network
        objs = sim.geometries.values()
    else:
        raise Exception('Setting indices does not apply to received object')

    # Ensure indices are not already assigned to another object
    temp = sp.zeros(boss._count(element=element), dtype=int)
    for item in objs:
        inds = boss._get_indices(element=element, labels=item.name)
        temp[inds] += 1
    temp[indices] += 1  # Increment proposed indices
    if sp.any(temp[indices] > 1):
        raise Exception('Some of the given '+element+' are already ' +
                        'assigned to an existing object')

    # Create new 'all' label for new size
    new_len = obj._count(element=element) + sp.size(indices)
    obj.update({element+'.all': sp.ones((new_len, ), dtype=bool)})

    # Set indices in boss dictionary
    inds_orig = boss._get_indices(element=element, labels=obj.name)
    if element+'.'+obj.name not in boss.keys():
        boss[element+'.'+obj.name] = False
    boss[element+'.'+obj.name][indices] = True
    inds_new = boss._get_indices(element=element, labels=obj.name)

    # Increase size of labels (add False at new indices)
    labels = obj.labels()
    labels.remove(element+'.all')
    for item in labels:
        if item.split('.')[0] == element:
            boss[element+'.'+'blank'] = False
            boss[element+'.'+'blank'][inds_orig] = obj[item]
            obj[item] = boss[element+'.'+'blank'][inds_new]
    boss.pop(element+'.'+'blank', None)


def drop_indices(obj, element, indices):
    sim = obj.simulation
    element = obj._parse_element(element, single=True)
    indices = obj._parse_indices(indices)

    # Adapt the method depending on type of object received
    if 'GenericPhysics' in obj.mro():
        boss = sim.find_phase(obj)
    elif 'GenericGeometry' in obj.mro():
        boss = sim.network
    else:
        raise Exception('Setting indices does not apply to received object')
    # Change the labeling in the boss object
    boss[element+'.'+obj.name][indices] = False
    # Convert network indices to obj-specific indices
    obj_inds = boss._map(element=element,
                         locations=indices,
                         target=obj)
    keep = ~obj._tomask(indices=obj_inds, element=element)
    for item in list(obj.keys()):
        if item.split('.')[0] == element:
            temp = obj[item][keep]
            obj.update({item: temp})
