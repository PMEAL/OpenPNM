'''
###############################################################################
Controller:  Overall simulation controller class
###############################################################################
'''
import pprint as pp

class Controller(dict):
    r"""

    """

    def __init__(self):
        pass

    def __str__(self):
        for item in self.keys():
            print("{a:<20s} {b:<35s}".format(a=self[item].__class__.__name__, b=item))
        return ''
        
    def network(self):
        return self._get_objects(obj_type='GenericNetwork')
        
    def geometries(self):
        return self._get_objects(obj_type='GenericGeometry')
        
    def phases(self):
        return self._get_objects(obj_type='GenericPhase')
        
    def physics(self):
        return self._get_objects(obj_type='GenericPhysics')

    def algorithms(self):
        return self._get_objects(obj_type='GenericAlgorithm')
        
    def _get_objects(self,obj_type):
        temp = []
        for obj in self.keys():
            mro = [item.__name__ for item in self[obj].__class__.__mro__]
            if obj_type in mro:
                temp.append(self[obj])
        return temp

    def drop(self,obj):
        r'''
        '''
        name = obj.name
        # Remove label arrays from all other objects
        for item in self.keys():
            self[item].pop('pore.'+name,None)
            self[item].pop('throat.'+name,None)
        # Set object's simulation to an empty dict
        self[name]._sim = {}
        # Remove object from simulation dict
        del self[name]
        
    def add(self,obj):
        r'''
        '''
        # Update object's simulation attribute
        obj.simulation = self
        # Add object to simulation dict
        self.update({obj.name: obj})

    @classmethod
    def _add_logger(cls,**kwargs):
        import logging as _logging
        _logging.basicConfig(level=_logging.ERROR,
                             format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                             datefmt='%m-%d %H:%M',
                             )
        if 'loggername' in kwargs.keys():
            cls._logger = _logging.getLogger(kwargs['loggername'])
        else:
            cls._logger = _logging.getLogger(cls.__class__.__name__)
        if 'loglevel' in kwargs.keys():
            cls._loglevel = kwargs['loglevel']


if __name__ == '__main__':
    sim = Controller()












