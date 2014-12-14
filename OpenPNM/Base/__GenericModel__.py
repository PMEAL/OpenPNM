'''
###############################################################################
GenericModel:  Abstract Class for Defining and Manipulating a Model
###############################################################################
'''
import scipy as sp
from OpenPNM.Base import logging, Controller
logger = logging.getLogger()

class GenericModel(dict):
    r"""
    Accepts a model from the OpenPNM model library, as well as all required and
    optional argumnents, then wraps it in a custom dictionary with various 
    methods for working with the models.

    """
    def __init__(self,**kwargs):
        self.update(**kwargs)

    def __call__(self):
        return self['model'](**self)
        
    def __str__(self):
        header = '-'*60
        print(header)
        print(self['model'].__module__+'.'+self['model'].__name__)
        print(header)
        print("{a:<20s} {b}".format(a='Argument Name',b='Value'))
        print(header)
        for item in self.keys():
            if item not in ['model','network','geometry','phase','physics','propname']:
                print("{a:<20s} {b}".format(a=item, b=self[item]))
        print(header)
        return ' '
        
    def regenerate(self):
        r'''
        
        '''
        self['func']()
        

if __name__ == '__main__':
    import OpenPNM.Geometry.models.pore_misc as mods
    a = GenericModel()



