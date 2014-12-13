'''
###############################################################################
GenericModel:  Abstract Class
###############################################################################
'''
import scipy as sp
from OpenPNM.Base import logging, Controller
logger = logging.getLogger()

class GenericModel(dict):
    r"""
    Accepts a model from the OpenPNM model library, as well as all required and
    optional argumnents, then wrapp

    """
    def __init__(self,model,**kwargs):
        self['func'] = model
        self.update(kwargs)

    def __call__(self):
        return self['func'](**self)
        
    def __repr__(self):
        f = self['func']
        header = '-'*60
        print(header)
        print(f.__module__+'.'+f.__name__)
        print(header)
        print("{a:<20s} {b}".format(a='Variable Name',b='Value'))
        print(header)
        for item in self.keys():
            if item not in ['network','geometry','phase','physics','propname']:
                print("{a:<20s} {b}".format(a=item, b=self[item]))
        print(header)
        return ''

if __name__ == '__main__':
    import OpenPNM.Geometry.models.pore_misc as mods
    a = GenericModel(model=mods.constant,value=3)



