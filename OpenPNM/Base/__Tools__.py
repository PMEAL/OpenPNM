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


class PrintableDict(_odict):
    def __repr__(self):
        text = dict(self).__str__()
        return text

    def __str__(self):
        horizontal_rule = '-' * 60
        lines = [horizontal_rule]
        lines.append('{0:<25s} {1}'.format('key', 'value'))
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
    def _get_health(self):
        health = True
        for item in list(self.keys()):
            if self[item] != []:
                health = False
        return health

    health = property(fget=_get_health)
