"""
###############################################################################
Tools:  Useful classes for use throughout the project
###############################################################################
"""
import scipy as _sp
from collections import OrderedDict as _odict

class PrintableList(list):
    def __str__(self):
        count = 0
        header = '-'*60
        print('\n')
        print(header)
        self.sort()
        for item in self:
            count = count + 1
            print(count,'\t: ',item)
        return header

class PrintableDict(_odict):
    def __repr__(self):
        text = dict(self).__str__()
        return text

    def __str__(self):
        header = '-'*60
        print(header)
        print("{a:<25s} {b:<25s}".format(a='key', b='value'))
        print(header)
        for item in self.keys():
            print("{a:<25s} {b}".format(a=item, b=self[item]))
        print(header)
        return ''
        
class AttributVeiew(object):
    def __init__(self, d):
        temp = {}
        for item in d:
            if type(d[item][0]) == _sp.bool_:
                key = 'label_'+item.replace('.','_')
            else:
                key = 'prop_'+item.replace('.','_')
            temp[key] =d[item]
        self.__dict__ = temp

class ClonedCore(dict):
    def __init__(self,obj):
        self.update(obj)
        self.name = obj.name

class HealthDict(PrintableDict):
    r"""
    This class adds a 'health' check to a standard dictionary.  This check 
    looks into the dict values, and considers empty lists as healthy and all 
    else as unhealthy.  If one or more entries is 'unhealthy' the health method
    returns False.
    """            
    def _get_health(self):
        health = True
        for item in self.keys():
            if self[item] != []:
                health = False
        return health
    
    health = property(fget=_get_health)
        
        
        
        
        
        