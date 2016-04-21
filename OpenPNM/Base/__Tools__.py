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
