'''
###############################################################################
Tools:  Useful classes for use throughout the project
###############################################################################
'''
import scipy as _sp

class PrintableList(list):
    def __str__(self):
        count = 0
        header = '-'*79
        print('\n')
        print(header)
        self.sort()
        for item in self:
            count = count + 1
            print(count,'\t: ',item)
        return header

class PrintableDict(dict):
    def __str__(self):
        import pprint
        header = '-'*79
        print('\n')
        print(header)
        pprint.pprint(self)
        print(header)
        return ''
        
class Attributeiew(object):
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
